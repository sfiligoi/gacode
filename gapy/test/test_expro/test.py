import gapy
import numpy as np

#gapy.expro.expro_read_legacy()
#gapy.expro.expro_write()

gapy.expro.expro_read('./')

nion = gapy.expro.expro_n_ion

print 'nexp    ',gapy.expro.expro_n_exp
print 'name    ',gapy.expro.expro_name[0]
print 'torfluxa',gapy.expro.expro_torfluxa
print 'ni      ',gapy.expro.expro_ni[0,:]

# Derived

print 'bunit   ',gapy.expro.expro_bunit
