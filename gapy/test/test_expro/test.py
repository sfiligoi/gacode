import gapy
import numpy as np

gapy.expro.expro_read('input.gacode')

nion = gapy.expro.expro_n_ion

print 'nexp    ',gapy.expro.expro_n_exp
print 'torfluxa',gapy.expro.expro_torfluxa
print 'ni      ',gapy.expro.expro_ni[0,:]

# Derived

print 'bunit   ',gapy.expro.expro_bunit
