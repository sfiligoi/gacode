from gacode import expro
import numpy as np
import sys

# Function to decode the insane string returned by gacode/f2py
def gapystr(s):
   if sys.version_info[0] == 2:
      return str(s).split()
   else:
      return str(s,'utf-8').split()

expro.expro_read('../data/input.gacode')

nion = expro.expro_n_ion

for i in range(nion):
    print(gapystr(expro.expro_name)[i]+gapystr(expro.expro_type)[i])

print('nexp  : '+str(expro.expro_n_exp))
print('bunit : '+str(expro.expro_bunit))
