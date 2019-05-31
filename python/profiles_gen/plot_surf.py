import gapy
import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text',usetex=True)
rc('font',size=18)

# Number of theta-points for plotting
narc = 256

surf  = sys.argv[1]
n     = int(sys.argv[2])
ftype = sys.argv[3]

# Read profiles
gapy.expro.expro_read('input.gacode')
nexp = int(gapy.expro.expro_n_exp) 
nfourier = int(gapy.expro.expro_nfourier)

print('nexp     = {:d}'.format(nexp))
print('nfourier = {:d}'.format(nfourier))

fig = plt.figure(figsize=(8,12))
ax = fig.add_subplot(111,aspect='equal')
ax.set_xlabel(r'$R$')
ax.set_ylabel(r'$Z$')

t = 2*np.pi*np.linspace(0,1,narc)

if n > 0:
   rlist = np.arange(0,nexp,n)
else:
   rlist = np.arange(abs(n),abs(n)+1)

geo_ar = gapy.expro.expro_geo[0,:,:]
geo_br = gapy.expro.expro_geo[1,:,:]
geo_az = gapy.expro.expro_geo[2,:,:]
geo_bz = gapy.expro.expro_geo[3,:,:]

if surf == 'msurf' or surf == 'surf':
   # Miller geometry flux-surfaces
   for i in rlist:
      bigr = gapy.expro.expro_rmaj[i]  
      bigz = gapy.expro.expro_zmag[i]   
      r    = gapy.expro.expro_rmin[i]
      d    = gapy.expro.expro_delta[i]
      k    = gapy.expro.expro_kappa[i]
      z    = gapy.expro.expro_zeta[i]
      
      x = bigr+r*np.cos(t+np.arcsin(d)*np.sin(t))
      y = bigz+k*r*np.sin(t+z*np.sin(2*t))

      if i == rlist[0]:
         ax.plot(x,y,'-k',linewidth=1,label=r'$\mathrm{Miller}$')
      else:
         ax.plot(x,y,'-k',linewidth=1)


if surf == 'fsurf' or surf == 'surf':
   # Fourier geometry flux-surfaces 
   for i in rlist:
      x = geo_ar[0,i]/2
      y = geo_az[0,i]/2
      for j in range(nfourier):
         p = j+1
         ar = geo_ar[p,i]
         br = geo_br[p,i]
         az = geo_az[p,i]
         bz = geo_bz[p,i]
         x = x+ar*np.cos(p*t)+br*np.sin(p*t)
         y = y+az*np.cos(p*t)+bz*np.sin(p*t)

      if i == rlist[0]:
         ax.plot(x,y,'-b',linewidth=1,label=r'$\mathrm{Fourier}~'+str(nfourier)+'$')
      else:
         ax.plot(x,y,'-b',linewidth=1)

# LCFS
i = nexp-1
x = geo_ar[0,i]/2
y = geo_az[0,i]/2

for j in range(nfourier):
   p = j+1
   ar = geo_ar[p,i]
   br = geo_br[p,i]
   az = geo_az[p,i]
   bz = geo_bz[p,i]
   x = x+ar*np.cos(p*t)+br*np.sin(p*t)
   y = y+az*np.cos(p*t)+bz*np.sin(p*t)

ax.plot(x,y,'-m',linewidth=2)

ax.legend()
if ftype == 'screen':
   plt.show()
else:
   outfile = ftype
   plt.savefig(outfile)
