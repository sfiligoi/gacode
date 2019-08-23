# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import gapy
import numpy as np

gapy.vpro.vpro_read_legacy()

gapy.vpro.vpro_write()
gapy.vpro.vpro_init(0)

gapy.vpro.vpro_read()

print('nexp    ',gapy.vpro.nexp)
print('bt_exp  ',gapy.vpro.bt_exp)
print('arho_exp',gapy.vpro.arho_exp)
print('ni      ',gapy.vpro.ni[:,0])

gapy.vpro.vpro_compute_derived()

print('bunit   ',gapy.vpro.bunit)
