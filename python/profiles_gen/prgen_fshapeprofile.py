import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from prgen_geqdsk import *
from prgen_contour import *
from prgen_fshape import *

if len(sys.argv) > 1:
   gfile = sys.argv[1]
   nrz   = int(sys.argv[2])
   npsi  = int(sys.argv[3])
   nf    = int(sys.argv[4])
else:
   print('Usage: python prgen_fshapeprofile.py <gfile> <nrz> <npsi> <nf>')
   sys.exit()

efit = prgen_geqdsk(gfile)
n_arc = 512

ri,zi,psi,q,p,fpol = prgen_contour(efit,nrz=nrz,levels=npsi,psinorm=0.999,narc=n_arc,quiet=False)

ari = np.zeros([nf+1,npsi])
bri = np.zeros([nf+1,npsi])
azi = np.zeros([nf+1,npsi])
bzi = np.zeros([nf+1,npsi])

for i in range(npsi-1):
   r=ri[:,i+1] ; z=zi[:,i+1]
   ar,br,az,bz = prgen_fshape(r,z,nf)
   ari[:,i+1] = ar[:]
   bri[:,i+1] = br[:]
   azi[:,i+1] = az[:]
   bzi[:,i+1] = bz[:]

u = ari
u = np.append(u,bri)
u = np.append(u,azi)
u = np.append(u,bzi)
u.tofile('fluxfit.geo')
