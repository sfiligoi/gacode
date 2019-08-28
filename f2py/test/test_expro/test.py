from gacode import expro
import numpy as np

expro.expro_read('../data/input.gacode')

nion = expro.expro_n_ion

print('nexp  : '+str(nion))
print('names : ',expro.expro_name)
print('bunit : ',expro.expro_bunit)
