import sys
from gacode import expro

# Function to decode the insane string returned by gacode/f2py
def gapystr(s):
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


class gapy(dict):
    def __init__(self, filename, input_profiles_compatibility_mode=True):
        self.filename = filename
        self.input_profiles_compatibility_mode = input_profiles_compatibility_mode
        self.load

    def load(self):
        expro.expro_read(self.filename)
        list = gapystr(expro.expro_list)

        # Define Gapy class members corresponding to list[]
        for item in list:
            self[item] = getattr(expro, 'expro_' + item)

        # Species name and type
        self['name'] = gapystr(expro.expro_name)
        self['type'] = gapystr(expro.expro_type)

        # refactor output to be compatible with input.profiles format
        if self.input_profiles_compatibility_mode:
            for item in ['n_exp', 'shot', 'time']:
                tmp = self[item]
                del self[item]
                self[item.upper()] = tmp
            self['Te'] = self['te']
            del self['te']
            self['Ti'] = self['ti']
            del self['ti']
            for quantity in ['ni', 'Ti', 'vpol', 'vtor']:
                if quantity in self:
                    for k in range(self['n_ion']):
                        self[quantity + '_%d' % (k + 1)] = self[quantity][k]
                    del self[quantity]
            self['IONS'] = {}
            for k in range(self['n_ion']):
                self['IONS'][k + 1] = [self['name'][k], self['z'][k], self['mass'][k], self['type'][k]]
            del self['name']
            del self['type']
            del self['z']
            del self['mass']
            del self['n_ion']
            for item in self:
                if isinstance(self[item], numpy.ndarray) and not len(self[item].shape):
                    self[item] = self[item].item()

    def sort(self):
        keys = sorted(list(self.keys()), key=lambda x: x.lower())
        for item in ['torfluxa', 'current', 'rcentr', 'bcentr', 'TIME', 'SHOT', 'N_EXP', 'IONS']:
            if item in keys:
                keys.pop(keys.index(item))
                keys = [item] + keys
        tmp = {}
        tmp.update(self)
        self.clear()
        for item in keys:
            self[item] = tmp[item]
