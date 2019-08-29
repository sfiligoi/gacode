from gacode import expro
import numpy
import sys

# Function to decode the insane string returned by gacode/f2py
def gapystr(s):
   if sys.version_info[0] == 2:
      return str(s).split()
   else:
      return str(s,'utf-8').split()

class Gapy(dict):

    def __init__(self, filename, input_profiles_compatibility_mode=True):
        expro.expro_read(filename)

        # input.gacode quantities
        list = ['n_exp','n_ion','mass','z','torfluxa','rvbv','ipa',
                'rho','rmin','polflux','q','w0','rmaj','zmag',
                'kappa','delta','zeta','ne','ni','te','ti','ptot',
                'johm','jbs','jrf','jnb','jbstor','sigmapar',
                'z_eff','vpol','vtor',
                'qohme','qbeame','qbeami','qrfe','qrfi','qfuse','qfusi',
                'qbrem','qsync','qline','qei','qione','qioni','qcxi','qpar','qmom']
        # Define Gapy class members corresponding to list[]
        for item in list:
           self[item] = getattr(expro,'expro_'+item)

        # Species name and type
        self['name'] = gapystr(expro.expro_name)
        self['type'] = gapystr(expro.expro_type)

        # Selected derived quantities
        list = ['bunit','gamma_e','gamma_p','s','drmaj','dzmag',
                'sdelta','skappa','szeta','dlnnedr','dlntedr','w0p',
                'vol','volp','cs','rhos','nuee']

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

