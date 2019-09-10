import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from prgen_geqdsk import *
from prgen_contour import *
from prgen_shape import *

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

pnorm = np.sqrt((psi[:]-psi[0])/(psi[-1]-psi[0]))

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
   ax.set_xlim([0.9,1])
#   ax.plot(pnorm,ci[i,:],'-m',linewidth=1,alpha=1)
   ax.plot(pnorm,np.gradient(ci[i,:],pnorm),'-m',linewidth=1,alpha=1)

label=['\kappa','\delta','\zeta','s_3','s_4']
for i in range(nf+1):
   ax = fig.add_subplot(2,nf+1,i+1+nf+1)
   ax.set_xlabel(r'$\psi$')
   ax.set_title(r'$'+label[i]+'$')
   ax.grid(which="both",ls=":")
   ax.set_xlim([0.9,1])
   if i > 0:
#      ax.plot(pnorm,si[i,:],'-m',linewidth=1,alpha=1)
      ax.plot(pnorm,np.gradient(si[i,:],pnorm),'-m',linewidth=1,alpha=1)
   else:
#      ax.plot(pnorm,xi[1,:],'-m',linewidth=1,alpha=1)
      ax.plot(pnorm,np.gradient(xi[1,:],pnorm),'-m',linewidth=1,alpha=1)
      
      
plt.tight_layout()
plt.show()
