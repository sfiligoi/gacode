import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from prgen_geqdsk import *
from prgen_contour import *
from prgen_shape import *

def extrap(x,u):
   m = (u[5]-u[4])/(x[5]-x[4])
   b = u[5]-m*x[5]
   u[0] = b
   u[1] = m*x[1]+b
   u[2] = m*x[2]+b
   u[3] = m*x[3]+b
   return u

def zero(x,u):
   r = u[4]/x[4]
   u[0] = 0.0
   u[1] = x[1]*r
   u[2] = x[2]*r
   u[3] = x[3]*r
   return u

def iring(x,u,xm):
   i = np.argmin(np.abs(x-xm))
   x0 = x[i]
   y0 = u[i]
   d = 0.0
   for j in range(8):
      d = y0*(1.0+d-x0)

   z = d/(1+d-x)
   return z,i

def ising(x,u,xm):
   i = np.argmin(np.abs(x-xm))
   x0 = x[i]
   y0 = u[i]
   r  = 1/y0
   d = 0.0
   for j in range(8):
      d = (1.0+d-x0)**r
      
   z = np.log(1+d-x)/np.log(d)
   return z,i

if len(sys.argv) > 1:
   gfile = sys.argv[1]
   npsi  = int(sys.argv[2])
   nrz   = int(sys.argv[3])
   ix    = int(sys.argv[4])
   nf    = int(sys.argv[5])
else:
   print('Usage: python prgen_shapeprofile.py <gfile> <npsi> <nrz> <ix> <nf>')
   sys.exit()

efit = prgen_geqdsk(gfile)
n_arc = 512

ri,zi,psi,q,p,fpol = prgen_contour(efit,nrz=nrz,levels=npsi,psinorm=0.999,narc=n_arc,quiet=False)

pnorm = ((psi[:]-psi[0])/(psi[-1]-psi[0]))
rnorm = np.sqrt(pnorm)

# ci -> cosine terms
# si -> sine terms
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
xi[1,:] = extrap(rnorm,xi[1,:]) # rmaj
xi[2,:] = extrap(rnorm,xi[2,:]) # kappa
xi[3,:] = extrap(rnorm,xi[3,:]) # zmaj

si[1,:] = zero(rnorm,si[1,:]) # delta
si[2,:] = zero(rnorm,si[2,:]) # zeta
si[3,:] = zero(rnorm,si[3,:]) # s3

ci[0,:] = extrap(rnorm,ci[0,:]) # tilt
ci[1,:] = zero(rnorm,ci[1,:]) # c1
ci[2,:] = zero(rnorm,ci[2,:]) # c2
ci[3,:] = zero(rnorm,ci[3,:]) # c3

# Repair near separatrix
si0 = np.zeros([nf+1,npsi])
ci0 = np.zeros([nf+1,npsi])
si0[:,:] = si[:,:]
ci0[:,:] = ci[:,:]

u = si0[1,:]
z,i = ising(rnorm,u/u[-1],0.9)
si[1,i:] = u[-1]*z[i:]

u = si0[2,:]
z,i = iring(pnorm,u/u[-1],0.81)
si[2,i:] = u[-1]*z[i:]

u = ci0[1,:]
z,i = iring(pnorm,u/u[-1],0.81)
ci[1,i:] = u[-1]*z[i:]

u = ci0[2,:]
z,i = iring(pnorm,u/u[-1],0.81)
ci[2,i:] = u[-1]*z[i:]

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
   ax.set_xlim([0.0,1])
   ax.plot(rnorm,xi[i,:],'-k',linewidth=1,alpha=1)

label=['c_0','c_1','c_2','c_3']
for i in range(4):
   ax = fig.add_subplot(3,4,i+5)
   ax.set_xlabel(r'$\psi$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0.0,1])
   u = ci[i,:] ; u0 = ci0[i,:]
   ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
   ax.plot(pnorm,u0,'-k',linewidth=1,alpha=1)

label=['-','\delta','-\zeta','s_3']
for i in range(4):
   ax = fig.add_subplot(3,4,i+9)
   ax.set_xlabel(r'$\psi$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0.0,1.0])
   if i > 0:
      u = si[i,:] ; u0 = si0[i,:]
      ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
      ax.plot(pnorm,u0,'-k',linewidth=1,alpha=1)
      
      
plt.tight_layout()
plt.show()
