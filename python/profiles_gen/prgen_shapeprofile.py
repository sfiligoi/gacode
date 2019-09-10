import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from prgen_geqdsk import *
from prgen_contour import *
from prgen_shape import *

def ring(x,u,xm):
   i = np.argmin(np.abs(x-xm))
   x0 = x[i]
   y0 = u[i]
   d = 0.0
   for j in range(8):
      d = y0*(1.0+d-x0)

   z = d/(1+d-x)
   return z

def iring(x,u,xm):
   i = np.argmin(np.abs(x-xm))
   x0 = x[i]
   y0 = u[i]
   d = 0.0
   for j in range(8):
      d = y0*(1.0+d-x0)

   z = d/(1+d-x)
   return z,i

def sing(x,u,xm):
   i = np.argmin(np.abs(x-xm))
   x0 = x[i]
   y0 = u[i]
   r  = 1/y0
   d = 0.0
   for j in range(8):
      d = (1.0+d-x0)**r
      
   z = np.log(1+d-x)/np.log(d)
   return z

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
#print(fpol)

pnorm = ((psi[:]-psi[0])/(psi[-1]-psi[0]))
rnorm = np.sqrt(pnorm)

# ci -> cosine terms
# si -> sine terms
# xi -> rmin,rmaj,kappa,delta

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
   ci[0,0] = 2*ci[0,1]-ci[0,2] # tilt
   xi[:,0] = 2*xi[:,1]-xi[:,2] # tilt
else:
   r=ri[:,ix] ; z=zi[:,ix]
   cr,sr,xr = prgen_shape(r,z,n_arc,nf,True)
   sys.exit()

# Repair functions
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

u = ci0[3,:]
z,i = iring(pnorm,u/u[-1],0.81)
ci[3,i:] = u[-1]*z[i:]

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

   
fig = plt.figure(figsize=(14,8))

label=['c_0','c_1','c_2','c_3','c_4']
for i in range(nf+1):
   ax = fig.add_subplot(2,nf+1,i+1)
   ax.set_xlabel(r'$\psi$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0.0,1])
   u = ci[i,:] ; u0 = ci0[i,:]
   ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
   ax.plot(pnorm,u0,'-k',linewidth=1,alpha=1)

label=['\kappa','\delta','\zeta','s_3','s_4']
for i in range(nf+1):
   ax = fig.add_subplot(2,nf+1,i+1+nf+1)
   ax.set_xlabel(r'$\psi$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0.0,1])
   if i > 0:
      u = si[i,:] ; u0 = si0[i,:]
      if i == 1:         
         ax.plot(rnorm,u,'-r',linewidth=1,alpha=1)
         ax.plot(rnorm,u0,'-k',linewidth=1,alpha=1)
      else:
         ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
         ax.plot(pnorm,u0,'-k',linewidth=1,alpha=1)
   else:
      ax.plot(pnorm,xi[1,:],'-k',linewidth=1)
      
      
plt.tight_layout()
plt.show()
