import sys
from gacode import expro
from gacodefuncs import *

class Gapy(dict):

    def __init__(self, filename, input_profiles_compatibility_mode=True):
        expro.expro_read(filename)
        list=gapystr(expro.expro_list)

        # Define Gapy class members corresponding to list[]
        for item in list:
           print('* '+item)
           self[item] = getattr(expro,'expro_'+item)

        # Species name and type
        self['name'] = gapystr(expro.expro_name)
        self['type'] = gapystr(expro.expro_type)

        if input_profiles_compatibility_mode:
            self['Te'] = self['te']
            del self['te']
            self['Ti'] = self['ti']
            del self['ti']
            for quantity in ['ni','Ti','vpol','vtor']:
                if quantity in self:
                    for k in range(self['n_ion']):
                        self[quantity + '_%d' % (k + 1)] = self[quantity][k]
                    del self[quantity]
            self['IONS'] = []
            for k in range(self['n_ion']):
                self['IONS'].append([self['name'][k], self['z'][k], self['mass'][k], self['type'][k]])
            del self['z']
            del self['mass']
            print(self['IONS'])

import cgyro
import gyro
import neo
import tgyro
from gacode import geo

