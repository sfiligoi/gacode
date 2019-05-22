import struct
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from gacodefuncs import *
from cgyro.data import cgyrodata
from mayavi import mlab 
try:
   import gapy
except:
   print('ERROR: (vis_torcut) Please build gapy/f2py library!')
   sys.exit()
   
ext      = sys.argv[1]
moment   = sys.argv[2]
species  = int(sys.argv[3])
nx       = int(sys.argv[4])
nz       = int(sys.argv[5])
istr     = sys.argv[6]
fmin     = sys.argv[7]
fmax     = sys.argv[8]
colormap = sys.argv[9]
font     = int(sys.argv[10])
legacy   = bool(int(sys.argv[11]))
dn       = int(sys.argv[12])
mag      = float(sys.argv[13])
nozonal  = bool(int(sys.argv[14]))

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

print('HINT: adjust -dn to match experimental dn (rho/a and Lx/a will shrink)')
print('Lx/rho = {:.2f}'.format(sim.length))
print('rho/a  = {:.4f}'.format(sim.rho/dn))
lovera = sim.length*sim.rho/dn*mag
print('Lx/a   = {:.4f}'.format(lovera))

ivec = time_vector(istr,nt)

s=ext.split('.')
if len(s) == 2:
   pre   = s[0]
   ftype = s[1]
else:
   pre = ''
   ftype = s[0]

#------------------------------------------------------------------------
# (r,theta)=(x,y) mesh setup 
#
if nth == 1:
   print('WARNING: (vis_torcut) Should use THETA_PLOT > 1 in CGYRO.')

if nx < 0:
   nx = 128
if nz < 0:
   nz = 128

x = np.zeros([nx])
z = np.zeros([nz])

for i in range(nx):
   x[i] = i*2*np.pi/(nx-1.0)/dn
for k in range(nz):
   z[k] = k*2*np.pi/(nz-1.0)-np.pi

xp = np.zeros([nx,nz])
yp = np.zeros([nx,nz])
zp = np.zeros([nx,nz])

for i in range(nx):
   r = sim.rmin+(dn*x[i]/(2*np.pi)-0.5)*lovera
   for k in range(nz):
      xp[i,k] = sim.rmaj+r*np.cos(z[k]+np.arcsin(sim.delta)*np.sin(z[k]))
      yp[i,k] = sim.zmag+sim.kappa*r*np.sin(z[k]+sim.zeta*np.sin(2*z[k]))
      zp[i,k] = 0.0

# Shape functions 
gapy.geo.signb_in=1 # fix
gapy.geo.geo_rmin_in=sim.rmin
gapy.geo.geo_rmaj_in=sim.rmaj
gapy.geo.geo_drmaj_in=sim.shift
gapy.geo.geo_zmag_in=sim.zmag
gapy.geo.geo_dzmag_in=sim.dzmag
gapy.geo.geo_q_in=sim.q
gapy.geo.geo_s_in=sim.shear
gapy.geo.geo_kappa_in=sim.kappa
gapy.geo.geo_delta_in=sim.delta
gapy.geo.geo_zeta_in=sim.zeta
gapy.geo.geo_s_kappa_in=sim.s_kappa
gapy.geo.geo_s_delta_in=sim.s_delta
gapy.geo.geo_s_zeta_in=sim.s_zeta
gapy.geo.geo_beta_star_in=sim.beta_star

gapy.geo.geo_interp(z,True)
if legacy:
   # s-alpha approximate (apparently used in legacy GYRO movies)
   # g1 -> q*theta
   # g2 -> theta 
   g1 = sim.q*z
   g2 = z
else:
   # Correct form of Clebsch angle expansion nu(r,theta) 
   g1 = -gapy.geo.geo_nu
   g2 = gapy.geo.geo_b*gapy.geo.geo_captheta/gapy.geo.geo_s_in/gapy.geo.geo_grad_r**2
   
if int(mag) == 0:
   showco=True
else:
   showco=False
if showco:
   # (Optional) plot of the geometry functions, then exit
   rc('text',usetex=True) ; rc('font',size=font)
   fig = plt.figure(figsize=(10,8))

   ax = fig.add_subplot(111)
   ax.plot(z/np.pi,g1/sim.q,label='g1')
   ax.plot(z/np.pi,g2,label='g2')
   ax.plot(z/np.pi,z,'--k')
   ax.set_xlim([-1,1])
   ax.legend()
   plt.show()
   sys.exit()
   
#------------------------------------------------------------------------

# Get filename and tags 
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists 
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print('INFO: (vis_torcut) Found binary data in '+fdata) 
else:
    print('ERROR: (vis_torcut) No data for -moment '+moment+' exists.  Try -moment phi')
    sys.exit()

if isfield:
    n_chunk = 2*nr*nn*nth
else:
    n_chunk = 2*nr*ns*nn*nth
   
# This is the logic to generate a frame
def frame():
   
   if isfield:
      a = np.reshape(aa,(2,nr,nth,nn),order='F')
   else:
      a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')
      
   mlab.figure(size=(900,900),bgcolor=(1,1,1))
   if isfield:
      c = a[0,:,:,:]+1j*a[1,:,:,:]
   else:
      c = a[0,:,:,species,:]+1j*a[1,:,:,species,:]
                
   if nozonal and nn > 1:
      c[:,:,0] = 0.0

   f = np.zeros([nx,nz],order='F')
   gapy.torcut(dn,sim.m_box,sim.q,sim.thetap,g1,g2,c,f)

   if fmin == 'auto':
      f0=np.min(f)
      f1=np.max(f)
   else:
      f0=float(fmin)
      f1=float(fmax)

   image = mlab.mesh(xp,yp,zp,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)
   # View from positive z-axis
   mlab.view(azimuth=0, elevation=0)
   print('INFO: (vis_torcut) min={:.3f} | max={:.3f}'.format(f0,f1))
   
   if ftype == 'screen':
      mlab.show()
   else:
      # Filename uses frame number 
      mlab.savefig(pre+str(i)+'.'+ftype)
      # Close each time to prevent memory accumulation
      mlab.close()
                

PREC='f' ; BIT=4

# Open binary file
fbin = open(fdata,'rb') 

i = 0
while True:
   try:
      aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
   except:
      sys.exit()

   i += 1
   print('INFO: (vis_torcut) Time index {:d} '.format(i))
   if i in ivec:
      frame()
   if i == max(ivec):
      sys.exit()
