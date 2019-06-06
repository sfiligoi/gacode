from pygacode import expro
import numpy as np

expro.expro_read('../data/input.gacode')

nion = expro.expro_n_ion

print('names : ',expro.expro_name)
print('bunit : ',expro.expro_bunit)
