from gacode import expro
import numpy as np
import sys

expro.expro_read('../data/input.gacode')

nion = expro.expro_n_ion

#for i in range(nion):
#    print(gapystr(expro.expro_name)[i]+gapystr(expro.expro_type)[i])

print('nexp  : '+str(expro.expro_n_exp))
print('bunit : '+str(expro.expro_bunit))
