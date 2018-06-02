import gapy
import numpy as np

gapy.vpro.vpro_read()

print 'nexp    ',gapy.vpro.nexp
print 'bt_exp  ',gapy.vpro.bt_exp
print 'arho_exp',gapy.vpro.arho_exp
print 'ni      ',gapy.vpro.ni[:,0]
