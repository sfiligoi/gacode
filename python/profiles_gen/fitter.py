import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from gacode_eqdsk import *
from gacode_polflux import *
from gacode_fit import *

if len(sys.argv) > 1:
   gfile = sys.argv[1]
   npsi  = int(sys.argv[2])
   nrz   = int(sys.argv[3])
   ix    = int(sys.argv[4])
   nf    = int(sys.argv[5])
else:
   print('Usage: python fitter.py <gfile> <npsi> <ix> <nf>')
   sys.exit()

EQDSK = geqdsk(gfile)
n_arc = 512

ri,zi,psi,q,p = polfluxcontour(EQDSK,nrz=nrz,levels=npsi,psinorm=0.9999,narc=n_arc,quiet=False)

pnorm = (psi[1:]-psi[0])/(psi[-1]-psi[0]) 

if ix < 1:
   ci = np.zeros([npsi,nf+1])
   si = np.zeros([npsi,nf+1])
   xi = np.zeros([npsi,4])
   for i in range(npsi-1):
      r=ri[:,i+1] ; z=zi[:,i+1]
      cr,sr,xr = fit(r,z,n_arc,nf,False)
      ci[i+1,:] = cr[:]
      si[i+1,:] = sr[:]
      xi[i+1,:] = xr[:]
   ci[0,0] = ci[1,0] # tilt
   xi[0,0] = 0.0 # rmin
   xi[0,1] = xi[1,1] # rmaj
   xi[0,2] = xi[1,2] # kappa
   xi[0,3] = xi[1,3] # zmaj
   
   
else:
   r=ri[:,ix] ; z=zi[:,ix]
   cr,sr = fit(r,z,n_arc,nf,True)
   sys.exit()

if ix == 0:
   f=open('out.dim','w')
   f.write(str(npsi))
   f.close()
   np.stack((si[:,0],
             si[:,1],
             si[:,2],
             ci[:,0],
             ci[:,1],
             ci[:,2],
             xi[:,0],
             xi[:,1],
             xi[:,2],
             xi[:,3],
             psi,
             q,
             p)).tofile('out.data')
   sys.exit()
   
#--------------------------------------------------------------------
# PLOTTING
#--------------------------------------------------------------------
# Latex fonts
rc('text',usetex=True)
rc('font',size=18)

fig = plt.figure(figsize=(12,8))

# tilt (c0)
ax = fig.add_subplot(231)
ax.set_xlabel(r'$\psi$')
ax.set_ylabel(r'$\mathrm{tilt}~(c_0)$')
ax.grid(which="both",ls=":")
ax.set_xlim([0,1])

# Data
ax.plot(pnorm,ci[:,0],'-k',linewidth=1,alpha=1)

# egg (c1)
ax = fig.add_subplot(232)
ax.set_xlabel(r'$\psi$')
ax.set_ylabel(r'$\mathrm{egg}~(c_1)$')
ax.grid(which="both",ls=":")
ax.set_xlim([0,1])

# Data
ax.plot(pnorm,ci[:,1],'-k',linewidth=1,alpha=1)

# egg (c2)
ax = fig.add_subplot(233)
ax.set_xlabel(r'$\psi$')
ax.set_ylabel(r'$\mathrm{shear}~(c_2)$')
ax.grid(which="both",ls=":")
ax.set_xlim([0,1])

# Data
ax.plot(pnorm,ci[:,2],'-k',linewidth=1,alpha=1)

# tilt (s0)
ax = fig.add_subplot(234)
ax.set_xlabel(r'$\psi$')
ax.set_ylabel(r'NONE')
ax.grid(which="both",ls=":")
ax.set_xlim([0,1])

# Data
ax.plot(pnorm,si[:,0],'-k',linewidth=1,alpha=1)

# egg (s1)
ax = fig.add_subplot(235)
ax.set_xlabel(r'$\psi$')
ax.set_ylabel(r'$\delta~(s_1)$')
ax.grid(which="both",ls=":")
ax.set_xlim([0,1])

# Data
ax.plot(pnorm,si[:,1],'-k',linewidth=1,alpha=1)

# egg (s2)
ax = fig.add_subplot(236)
ax.set_xlabel(r'$\psi$')
ax.set_ylabel(r'$\zeta~(s_2)$')
ax.grid(which="both",ls=":")
ax.set_xlim([0,1])

# Data
ax.plot(pnorm,si[:,2],'-k',linewidth=1,alpha=1)

plt.tight_layout()
plt.show()
