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
   npsi  = int(sys.argv[2])
   ix   = int(sys.argv[3])
   nrz   = int(sys.argv[4])
   nf    = int(sys.argv[5])
else:
   print('Usage: python prgen_fshapeprofile.py <gfile> <npsi> <ix> <nrz> <nf>')
   sys.exit()

#efit = prgen_geqdsk(gfile)
n_arc = 512

#ri,zi,psi,q,p,fpol = prgen_contour(efit,nrz=nrz,levels=npsi,psinorm=0.999,narc=n_arc,quiet=False)

#r=ri[:,ix] ; z=zi[:,ix]

t = np.linspace(0,1,n_arc)*2*np.pi
r = 3.0+0.5*np.cos(t+0.4*np.sin(t))*(1.0+0.2*np.sin(t))
z = 0.05+1.5*0.5*np.sin(t)

ar,br,az,bz = prgen_fshape(r,z,nf)

fig = plt.figure(figsize=(18,9))

# PLOT contour
ax = fig.add_subplot(111,aspect='equal')
ax.set_xlabel(r'$R$')
ax.set_ylabel(r'$Z$')
ax.grid(which="both",ls=":")

# Data
ax.plot(r,z,'-k',linewidth=2,alpha=0.3)
# Parameterized contour
rp = 0.5*ar[0]
zp = 0.5*az[0]
for n in range(1,nf+1):
   rp = rp+ar[n]*np.cos(n*t)+br[n]*np.sin(n*t)
   zp = zp+az[n]*np.cos(n*t)+bz[n]*np.sin(n*t)
   
ax.plot(rp,zp,'-m',linewidth=1)

plt.show()


