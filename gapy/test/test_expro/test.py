import gapy
import numpy as np

gapy.vpro.vpro_read_legacy()
gapy.vpro.vpro_write()

gapy.vpro.vpro_read()

print 'nexp    ',gapy.vpro.expro_n_exp
print 'bt_exp  ',gapy.vpro.expro_b_ref
print 'arho_exp',gapy.vpro.expro_arho
print 'ni      ',gapy.vpro.expro_ni[0,:]

# Derived

print 'bunit   ',gapy.vpro.expro_bunit
