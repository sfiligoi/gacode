import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from prgen_geqdsk import *
from prgen_contour import *
from prgen_shape import *

repair = False

if len(sys.argv) > 1:
   gfile   = sys.argv[1]
   nrz     = int(sys.argv[2])
   npsi    = int(sys.argv[3])
   ix      = int(sys.argv[4])
   nfourier = int(sys.argv[5])
else:
   print('Usage: python prgen_shapeprofile.py <gfile> <nrz> <npsi> <ix> <nfourier>')
   sys.exit()

efit  = prgen_geqdsk(gfile)
n_arc = 700
nf = 3

ri,zi,psi,q,p,fpol = prgen_contour(efit,nrz=nrz,levels=npsi,psinorm=0.995,narc=n_arc,quiet=False)

pnorm = ((psi[:]-psi[0])/(psi[-1]-psi[0]))
rnorm = np.sqrt(pnorm)

if 0==1:
   fig = plt.figure(figsize=(14,12))
   ax = fig.add_subplot(111)
   rd = np.gradient(ri[:,npsi//2])
   zd = np.gradient(zi[:,npsi//2])
   s = np.argmax(rd)

   # Shift elements so that first index at max(R).
   rd[0:-1] = np.roll(rd[0:-1],-s) ; rd[-1] = rd[0] 
   zd[0:-1] = np.roll(zd[0:-1],-s) ; zd[-1] = zd[0]
   ax.plot(rd-np.mean(rd))
   ax.plot(zd)
   plt.show()

if nfourier > 0:
   oldfourier(ri,zi,nfourier,rnorm)

# ci -> cos terms
# si -> sin terms
# xi -> rmin,rmaj,kappa,zmaj

if ix < 1:
   ci = np.zeros([nf+1,npsi])
   si = np.zeros([nf+1,npsi])
   xi = np.zeros([4,npsi])
   for i in range(npsi-1):
      r=ri[:,i+1] ; z=zi[:,i+1]
      cr,sr,xr = prgen_shape(r,z,n_arc,nf,False)
      ci[:,i+1] = cr[:]
      si[:,i+1] = sr[:]
      xi[:,i+1] = xr[:]
else:
   r=ri[:,ix] ; z=zi[:,ix]
   cr,sr,xr = prgen_shape(r,z,n_arc,nf,True)
   sys.exit()
   
# Repair functions near origin

xi[0,:] = zero(rnorm,xi[0,:]) # rmin
xi[1,:] = extrap(pnorm,xi[1,:]) # rmaj
xi[2,:] = extrap(rnorm,xi[2,:]) # kappa
xi[3,:] = extrap(rnorm,xi[3,:]) # zmaj

si[1,:] = zero(rnorm,si[1,:]) # delta
si[2,:] = zero(pnorm,si[2,:]) # zeta
si[3,:] = zero(pnorm,si[3,:]) # s3

ci[0,:] = extrap(pnorm,ci[0,:]) # tilt
ci[1,:] = zero(rnorm,ci[1,:]) # c1
ci[2,:] = zero(pnorm,ci[2,:]) # c2
ci[3,:] = zero(pnorm,ci[3,:]) # c3

# Repair near separatrix
si0 = np.zeros([nf+1,npsi]) ; si0[:,:] = si[:,:]
ci0 = np.zeros([nf+1,npsi]) ; ci0[:,:] = ci[:,:]
xi0 = np.zeros([4,npsi])    ; xi0[:,:] = xi[:,:]

if repair:
   u = si0[1,:] ; z,i = ising(rnorm,u/u[-1],0.92)
   si[1,i:] = u[-1]*z[i:]

   u = si0[2,:] ; z,i = iring(pnorm,u/u[-1],0.84)
   si[2,i:] = u[-1]*z[i:]

   u = ci0[0,:] ; z,i = iring(pnorm,(u-u[0])/(u[-1]-u[0]),0.84)
   ci[0,i:] = (u[-1]-u[0])*z[i:]+u[0]

   u = ci0[1,:] ; z,i = iring(pnorm,u/u[-1],0.84)
   ci[1,i:] = u[-1]*z[i:]

   u = ci0[2,:] ; z,i = iring(pnorm,u/u[-1],0.84)
   ci[2,i:] = u[-1]*z[i:]

   u = xi0[2,:] ; z,i = ising(pnorm,(u-u[0])/(u[-1]-u[0]),0.84)
   xi[2,i:] = (u[-1]-u[0])*z[i:]+u[0]

if ix == 0:
   f=open('out.dim','w')
   f.write(str(npsi)+'\n')
   f.write(str(nf+1)+'\n')
   f.write(str(efit['RCENTR'])+'\n')
   f.write(str(efit['BCENTR'])+'\n')
   f.write(str(efit['CURRENT']*1e-6)+'\n')
   f.close()
   u = psi
   u = np.append(u,q)
   u = np.append(u,p)
   u = np.append(u,si[:,:])
   u = np.append(u,ci[:,:])
   u = np.append(u,xi[:,:])
   u.tofile('out.data')
   sys.exit()

#--------------------------------------------------------------------
# PLOTTING
#--------------------------------------------------------------------
# Latex fonts
rc('text',usetex=True)
rc('font',size=18)

fig = plt.figure(figsize=(14,12))

label=['r','R','\kappa','Z']
for i in range(4):
   ax = fig.add_subplot(3,4,i+1)
   ax.set_xlabel(r'$\psi^{1/2}$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0,1])
   u = xi[i,:] ; u0 = xi0[i,:]
   ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
   ax.plot(pnorm,u0,'-k',linewidth=1,alpha=1)

label=['c_0','c_1','c_2','c_3']
for i in range(4):
   ax = fig.add_subplot(3,4,i+5)
   ax.set_xlabel(r'$\psi$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0,1])
   u = ci[i,:] ; u0 = ci0[i,:]
   ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
   ax.plot(pnorm,u0,'-k',linewidth=1,alpha=1)

label=['-','\delta','-\zeta','s_3']
for i in range(4):
   ax = fig.add_subplot(3,4,i+9)
   ax.set_xlabel(r'$\psi$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0,1])
   if i > 0:
      u = si[i,:] ; u0 = si0[i,:]
      ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
      ax.plot(pnorm,u0,'-k',linewidth=1,alpha=1)
            
plt.tight_layout()
plt.show()
