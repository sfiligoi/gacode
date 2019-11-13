import sys
from gacode import expro
import numpy
import re


def gapystr_get(s):
    if sys.version_info[0] == 2:
        u = []
        for i in range(len(s)):
            u.append(str(s[i]).strip())
    else:
        u = str(s, 'utf-8').split()
    out_u = []
    for s in u:
        if '\x00' in s:
            break
        else:
            out_u.append(s)
    return out_u


def gapystr_set(s):
    return [item.ljust(10) for item in s + [''] * 20][:20]


def sort_key(x):
    if re.match('[a-zA-Z]+_[0-9]+', x):
        x = x.split('_')
        x = '%s_%02d' % (x[0], int(x[1]))
    return x.lower()


class gapy(dict):
    def __init__(self, filename, input_profiles_compatibility_mode=True):
        self.filename = filename
        self.input_profiles_compatibility_mode = input_profiles_compatibility_mode
        self.load()

    def load(self):
        getattr(expro, 'expro_name')[:] = gapystr_set([''])
        getattr(expro, 'expro_type')[:] = gapystr_set([''])

        expro.expro_read(self.filename)
        lst = gapystr_get(expro.expro_list)

        # Define Gapy class members corresponding to list[]
        for item in lst:
            self[item] = getattr(expro, 'expro_' + item).copy()

        # Species name and type
        self['name'] = gapystr_get(expro.expro_name)
        self['type'] = gapystr_get(expro.expro_type)

        # refactor output to be compatible with input.profiles format
        if self.input_profiles_compatibility_mode:
            # capitalize some entries
            for item in ['n_exp', 'shot', 'time']:
                tmp = self[item]
                del self[item]
                self[item.upper()] = tmp
            # input.profiles Te and Ti are capitalized
            self['Te'] = self['te']
            del self['te']
            self['Ti'] = self['ti']
            del self['ti']
            self['omega0'] = self['w0']
            del self['w0']
            # input.profiles entries are split per ion species
            for item in list(self.keys()):
                if isinstance(self[item], numpy.ndarray) and len(self[item].shape) > 1:
                    for k in range(self['n_ion']):
                        self[item + '_%d' % (k + 1)] = self[item][k]
                    del self[item]
            # input.profiles ions infos are collected under the IONS structure
            self['IONS'] = {}
            for k in range(self['n_ion']):
                self['IONS'][k + 1] = [self['name'][k], self['z'][k], self['mass'][k], self['type'][k]]
            del self['name']
            del self['type']
            del self['z']
            del self['mass']
            del self['n_ion']
            # use floats and ints rather than 0D numpy arrays
            for item in self:
                if isinstance(self[item], numpy.ndarray) and not len(self[item].shape):
                    self[item] = self[item].item()

    def save(self):
        # write single arrays and collate items
        collated = {}
        for item in sorted(list(self.keys()), key=sort_key):
            if item == 'IONS':
                continue
            expro_item = item
            expro_value = self[item]
            if item in ['N_EXP', 'SHOT', 'TIME', 'Te']:
                expro_item = item.lower()
            if self.input_profiles_compatibility_mode and re.match('[a-zA-Z]+_[0-9]+', item):
                expro_item = item.split('_')
                collated.setdefault(expro_item[0].lower(), []).append(self[item])
                continue
            elif expro_item in ['name', 'type']:
                getattr(expro, 'expro_' + expro_item)[:] = gapystr_set(expro_value)
                continue
            if expro_item in expro.expro_list:
                setattr(expro, 'expro_' + expro_item, expro_value)
        # manage collated items
        for expro_item in collated:
            if expro_item in expro.expro_list:
                setattr(expro, 'expro_' + expro_item, collated[expro_item])
        # manage IONS
        if 'IONS' in self:
            for item_n, expro_item in enumerate(['name', 'z', 'mass', 'type']):
                expro_value = [self['IONS'][k][item_n] for k in self['IONS']]
                if expro_item in ['name', 'type']:
                    getattr(expro, 'expro_' + expro_item)[:] = gapystr_set(expro_value)
                    continue
                setattr(expro, 'expro_' + expro_item, expro_value)
            setattr(expro, 'expro_n_ion', len(self['IONS']))
        # write input.gacode file
        expro.expro_write(self.filename)

    def sort(self):
        # sort alphabetically independent of capitalization and consistent with subscript numbering
        keys = sorted(list(self.keys()), key=sort_key)
        # place some elements at the top of the input.profiles
        for item in ['torfluxa', 'current', 'rcentr', 'bcentr', 'TIME', 'SHOT', 'N_EXP', 'IONS']:
            if item in keys:
                keys.pop(keys.index(item))
                keys = [item] + keys
        # recreate dictionary with the correct sorting
        tmp = {}
        tmp.update(self)
        self.clear()
        for item in keys:
            self[item] = tmp[item]
