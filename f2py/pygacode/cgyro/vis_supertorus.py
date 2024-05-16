import struct
import sys
import time
import numpy as np
import os
from matplotlib import cm
from ..gacodefuncs import *
from .data import cgyrodata
from mayavi import mlab

try:
   from pygacode import geo
   from vis.vis import vis
   print('INFO: (vis_supertorus) Successfully imported vis')
except:
   print('ERROR: (vis_supertorus) Please type make in f2py/vis')
   sys.exit()

ext      = sys.argv[1]
moment   = sys.argv[2]
species  = int(sys.argv[3])
nx       = int(sys.argv[4])
nz       = int(sys.argv[5])
nphi     = int(sys.argv[6])
phimax   = float(sys.argv[7])
istr     = sys.argv[8]
fmin_in  = sys.argv[9]
fmax_in  = sys.argv[10]
colormap = sys.argv[11]
font     = int(sys.argv[12])
legacy   = bool(int(sys.argv[13]))
dn       = int(sys.argv[14])
mag      = float(sys.argv[15])
nozonal  = bool(int(sys.argv[16]))
onlyzonal = bool(int(sys.argv[17]))
px        = int(sys.argv[18])
py        = int(sys.argv[19])

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

lovera = sim.length*sim.rho/dn*mag

print('HINT: adjust -dn to match experimental dn (rho/a and Lx/a will shrink)')
print('Lx/rho = {:.2f}'.format(sim.length))
print('rho/a  = {:.4f}'.format(sim.rho/dn))
print('Lx/a   = {:.4f}'.format(lovera))
print('px = {}  py = {}'.format(px,py))
print('nx = {}  nz = {}  nphi = {}'.format(nx,nz,nphi))

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
start = time.time()

if nth == 1:
   print('WARNING: (vis_supertorus) Should use THETA_PLOT > 1 in CGYRO.')

if nx < 0:
   nx = 128
if nz < 0:
   nz = 128

x = np.linspace(0,2*np.pi,nx)/dn
z = np.linspace(0,2*np.pi,nz)-np.pi
   
xp = np.zeros([nx,nz])
yp = np.zeros([nx,nz])
zp = np.zeros([nx,nz])

# minor radius
r = sim.rmin+(dn*x/(2*np.pi)-0.5)*lovera

# MXH angle
a = z + (sim.shape_cos[0]
         + sim.shape_cos[1]*np.cos(z) 
         + sim.shape_cos[2]*np.cos(2*z) 
         + sim.shape_cos[3]*np.cos(3*z)
         + np.arcsin(sim.delta)*np.sin(z) 
         - sim.zeta*np.sin(2*z) 
         + sim.shape_sin[3]*np.sin(3*z))

# MESH
xp[:,:] = sim.rmaj+r[:,None]*np.cos(a[None,:])
yp[:,:] = sim.zmag+sim.kappa*r[:,None]*np.sin(z[None,:])
zp[:,:] = 0.0

print('MESH TIME = '+'{:.3e}'.format(time.time()-start)+' s.')

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
geo.geo_shape_cos0_in = sim.shape_cos[0]
geo.geo_shape_cos1_in = sim.shape_cos[1]
geo.geo_shape_cos2_in = sim.shape_cos[2]
geo.geo_shape_cos3_in = sim.shape_cos[3]
geo.geo_shape_cos4_in = 0.0
geo.geo_shape_cos5_in = 0.0
geo.geo_shape_cos6_in = 0.0
geo.geo_shape_sin3_in = sim.shape_sin[3]
geo.geo_shape_sin4_in = 0.0
geo.geo_shape_sin5_in = 0.0
geo.geo_shape_sin6_in = 0.0
geo.geo_shape_s_cos0_in = sim.shape_s_cos[0]
geo.geo_shape_s_cos1_in = sim.shape_s_cos[1]
geo.geo_shape_s_cos2_in = sim.shape_s_cos[2]
geo.geo_shape_s_cos3_in = sim.shape_s_cos[3]
geo.geo_shape_s_cos4_in = 0.0
geo.geo_shape_s_cos5_in = 0.0
geo.geo_shape_s_cos6_in = 0.0
geo.geo_shape_s_sin3_in = sim.shape_s_sin[3]
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

   mlab.figure(bgcolor=(1,1,1))
   if isfield:
      c = a[0,:,:,:]+1j*a[1,:,:,:]
   else:
      c = a[0,:,:,species,:]+1j*a[1,:,:,species,:]
                
   if nozonal and nn > 1:
      c[:,:,0] = 0.0
   if onlyzonal and nn > 1:
      c[:,:,1:] = 0.0
   
# This is the logic to generate an (r,theta) frame
def subframe_torcut(iframe):

   global c,fmin,fmax
   
   phi = -2*np.pi*phimax*iframe
   xpp = xp*np.cos(phi)
   ypp = yp
   zpp = xp*np.sin(phi)

   f = np.zeros([nx,nz],order='F')
   
   vis.torcut(dn,sim.m_box,sim.q,sim.thetap,g1-phi/dn,g2,c,f)

   if iframe == 0:
      if fmin_in == 'auto':
         fmin = np.min(f)
         fmax = np.max(f)
      else:
         fmin = float(fmin_in)
         fmax = float(fmax_in)
      print('INFO: (vis_supertorus) min={:.3f} | max={:.3f}'.format(fmin,fmax))
         
   image = mlab.mesh(xpp,ypp,zpp,scalars=f,
                     colormap=colormap,vmin=fmin,vmax=fmax,opacity=1.0)

# This is the logic to generate an (theta,phi) frame
def subframe_torside():

   global c,fmin,fmax
   
   phi = -np.linspace(0,2*np.pi*phimax,nphi)
   xpp = np.zeros([nz,nphi])
   ypp = np.zeros([nz,nphi])
   zpp = np.zeros([nz,nphi])
   
   f = np.zeros([2,nz,nphi],order='F')
   
   vis.torside(dn,sim.m_box,sim.q,phimax,sim.thetap,g1,g2,c,f)

   for i in [0,1]:
      xpp[:,:] = xp[-i,:,None]*np.cos(phi[None,:])
      ypp[:,:] = yp[-i,:,None]
      zpp[:,:] = xp[-i,:,None]*np.sin(phi[None,:])

      image = mlab.mesh(xpp[:,:],ypp[:,:],zpp[:,:],scalars=f[i,:,:],
                        colormap=colormap,vmin=fmin,vmax=fmax,opacity=1.0) 


PREC='f' ; BIT=4

print('ivec {}'.format(ivec))

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
         start = time.time()
         # cap at phi=0
         subframe_torcut(0)
         # cap at phi=-phimax
         subframe_torcut(1)
         # inner and outer body of torus
         subframe_torside()
         # View from positive z-axis
         mlab.view(azimuth=0,elevation=0)
         print('MAP TIME = '+'{:.3e}'.format(time.time()-start)+' s.')
         if ftype == 'screen':
            mlab.show()
         else:
            # Filename uses frame number 
            mlab.savefig(pre+str(i)+'.'+ftype,size=(px,py))
            # Close each time to prevent memory accumulation
            mlab.close()
         
      if i == max(ivec):
         sys.exit()
