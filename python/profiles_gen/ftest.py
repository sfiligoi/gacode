import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from prgen_geqdsk import *
from prgen_contour import *
from prgen_shape import *

t = np.linspace(0,1,101)*2*np.pi

r = 3+0.5*np.cos(t)+0.1*np.sin(t)
z = np.sin(t)
nf = 24

ar,br,az,bz = prgen_fshape(r,z,nf)

print(ar)
print(br)
print(az)
print(bz)
