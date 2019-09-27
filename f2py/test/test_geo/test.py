# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import gapy
import numpy as np

np.set_printoptions(precision=16)

theta=np.array([0.1,0.2])

gapy.geo.geo_interp(theta,True)

print(gapy.geo.geo_b)







