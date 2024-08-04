import struct
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from ..gacodefuncs import *
from .data import cgyrodata
from mayavi import mlab
try:
   from pygacode import vis
   from pygacode import geo
   print("INFO: (vis_supertorus) Successfully imported vis")
except:
   print("ERROR: (vis_supertorus) Please type 'pip install pygacode'")
   sys.exit()

ext      = sys.argv[1]
moment   = sys.argv[2]
species  = int(sys.argv[3])
nx       = int(sys.argv[4])
nz       = int(sys.argv[5])
nphi     = int(sys.argv[6])
phimax   = float(sys.argv[7])
istr     = sys.argv[8]
fmin     = sys.argv[9]
fmax     = sys.argv[10]
colormap = sys.argv[11]
font     = int(sys.argv[12])
legacy   = bool(int(sys.argv[13]))
dn       = int(sys.argv[14])
mag      = float(sys.argv[15])
nozonal  = bool(int(sys.argv[16]))
onlyzonal = bool(int(sys.argv[17]))

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
   print('WARNING: (vis_supertorus) Should use THETA_PLOT > 1 in CGYRO.')

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
      a = z[k] + (sim.shape_cos0
         + sim.shape_cos1*np.cos(z[k]) 
         + sim.shape_cos2*np.cos(2*z[k]) 
         + sim.shape_cos3*np.cos(3*z[k])
         + np.arcsin(sim.delta)*np.sin(z[k]) 
         - sim.zeta*np.sin(2*z[k]) 
         + sim.shape_sin3*np.sin(3*z[k]))
      xp[i,k] = sim.rmaj+r*np.cos(a)
      yp[i,k] = sim.zmag+sim.kappa*r*np.sin(z[k])
      zp[i,k] = 0.0

# Shape functions 
geo.signb_in=1 # fix
geo.geo_rmin_in=sim.rmin
geo.geo_rmaj_in=sim.rmaj
geo.geo_drmaj_in=sim.shift
geo.geo_zmag_in=sim.zmag
geo.geo_dzmag_in=sim.dzmag
geo.geo_q_in=sim.q
geo.geo_s_in=sim.shear
geo.geo_kappa_in=sim.kappa
geo.geo_delta_in=sim.delta
geo.geo_zeta_in=sim.zeta
geo.geo_s_kappa_in=sim.s_kappa
geo.geo_s_delta_in=sim.s_delta
geo.geo_s_zeta_in=sim.s_zeta
geo.geo_beta_star_in=sim.beta_star
geo.geo_shape_cos0_in = sim.shape_cos0
geo.geo_shape_cos1_in = sim.shape_cos1
geo.geo_shape_cos2_in = sim.shape_cos2
geo.geo_shape_cos3_in = sim.shape_cos3
geo.geo_shape_cos4_in = 0.0
geo.geo_shape_cos5_in = 0.0
geo.geo_shape_cos6_in = 0.0
geo.geo_shape_sin3_in = sim.shape_sin3
geo.geo_shape_sin4_in = 0.0
geo.geo_shape_sin5_in = 0.0
geo.geo_shape_sin6_in = 0.0
geo.geo_shape_s_cos0_in = sim.shape_s_cos0
geo.geo_shape_s_cos1_in = sim.shape_s_cos1
geo.geo_shape_s_cos2_in = sim.shape_s_cos2
geo.geo_shape_s_cos3_in = sim.shape_s_cos3
geo.geo_shape_s_cos4_in = 0.0
geo.geo_shape_s_cos5_in = 0.0
geo.geo_shape_s_cos6_in = 0.0
geo.geo_shape_s_sin3_in = sim.shape_s_sin3
geo.geo_shape_s_sin4_in = 0.0
geo.geo_shape_s_sin5_in = 0.0
geo.geo_shape_s_sin6_in = 0.0

geo.geo_interp(z,True)
if legacy:
   # s-alpha approximate (apparently used in legacy GYRO movies)
   # g1 -> q*theta
   # g2 -> theta 
   g1 = sim.q*z 
   g2 = z
else:
   # Correct form of Clebsch angle expansion nu(r,theta)
   g1 = -geo.geo_nu 
   g2 = geo.geo_b*geo.geo_captheta/geo.geo_s_in/geo.geo_grad_r**2
   
#------------------------------------------------------------------------

# Get filename and tags 
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists 
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print('INFO: (vis_supertorus) Found binary data in '+fdata) 
else:
    print('ERROR: (vis_supertorus) No data for -moment '+moment+' exists.  Try -moment phi')
    sys.exit()

if isfield:
    n_chunk = 2*nr*nn*nth
else:
    n_chunk = 2*nr*ns*nn*nth
    
def frame():
   global c
   
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
   if onlyzonal and nn > 1:
      c[:,:,1:] = 0.0
   
# This is the logic to generate an (r,theta) frame
def subframe_rt(iframe):
   global xp,yp,zp,g1,g2,c,fmin_all,fmax_all
   
   f = np.zeros([nx,nz],order='F')
   
   vis.torcut(dn,sim.m_box,sim.q,sim.thetap,g1,g2,c,f)

   if iframe == 0:
      if fmin == 'auto':
         fmin_all=np.min(f)
         fmax_all=np.max(f)
      else:
         fmin_all=float(fmin)
         fmax_all=float(fmax)
         
   image = mlab.mesh(xp,yp,zp,scalars=f,colormap=colormap,vmin=fmin_all,vmax=fmax_all,opacity=1.0)
   print('INFO: (vis_supertorus) min={:.3f} | max={:.3f}'.format(fmin_all,fmax_all))

# This is the logic to generate an (theta,phi) frame
def subframe_tp():
   global xp,yp,zp,g1,g2,c,g10,xp0,fmin_all,fmax_all

   dphi = 2*np.pi*phimax/(nphi-1)
   
   xpp    = np.zeros([2,nz,nphi])
   zpp    = np.zeros([2,nz,nphi])
   ypp    = np.zeros([2,nz,nphi])
   f_all  = np.zeros([2,nz,nphi],order='F')
   f      = np.zeros([2,nz],order='F')
   
   for j in range(nphi):
      #
      for k in range(nz):
         phi      = -j*dphi
         g1[k]    = g10[k] - phi
         xpp[0,k,j] = xp0[0,k] * np.cos(phi)
         zpp[0,k,j] = xp0[0,k] * np.sin(phi)
         ypp[0,k,j] = yp[0,k]
         xpp[1,k,j] = xp0[nx-1,k] * np.cos(phi)
         zpp[1,k,j] = xp0[nx-1,k] * np.sin(phi)
         ypp[1,k,j] = yp[nx-1,k]
      #
      f[:,:] = 0.0
      vis.torcut(dn,sim.m_box,sim.q,sim.thetap,g1,g2,c,f)
      for k in range(nz):
         f_all[0,k,j] = f[0,k]
         f_all[1,k,j] = f[1,k]

   for i in range(2):
      image = mlab.mesh(xpp[i,:,:],ypp[i,:,:],zpp[i,:,:],scalars=f_all[i,:,:],colormap=colormap,vmin=fmin_all,vmax=fmax_all,opacity=1.0) 

PREC='f' ; BIT=4

print("ivec")
print(ivec)

# Open binary file
with open(fdata,'rb') as fbin:
   i = 0
   while True:
      try:
         aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
      except:
         sys.exit()

      i += 1
      print('INFO: (vis_supertorus) Time index {:d} '.format(i))
      if i in ivec:
         frame()
         g10 = g1
         xp0 = xp
         ##################
         # cap at phi=0
         phi = 2*np.pi*0.0
         g1 = g10 - phi
         xp = xp0 * np.cos(phi)
         zp = xp0 * np.sin(phi)
         subframe_rt(0)
         ##################
         # cap at phi=phimax
         phi = -2*np.pi*phimax
         g1 = g10 - phi
         xp = xp0 * np.cos(phi)
         zp = xp0 * np.sin(phi)
         subframe_rt(1)
         ##################
         # inner and outer body of torus
         subframe_tp()
         ##################
         # View from positive z-axis
         mlab.view(azimuth=0, elevation=0)
         if ftype == 'screen':
            mlab.show()
         else:
            # Filename uses frame number 
            mlab.savefig(pre+str(i)+'.'+ftype)
            # Close each time to prevent memory accumulation
            mlab.close()
         
      if i == max(ivec):
         sys.exit()
