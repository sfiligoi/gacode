import struct
import sys
import numpy as np
import os
from matplotlib import rc
import matplotlib.pyplot as plt
from gacodefuncs import *
from cgyro.data import cgyrodata
from mayavi import mlab 
try:
   import gapy
   hasgapy = True
except:
   print ' BAD: (plot_viz) Please build gapy.so library!'
   sys.exit()
   
PREC='f' ; BIT=4

font=16
ftype = 'screen'
moment = 'phi'
species = 0
nx = 128
ny = 128
istr = "39"
fmin = 'auto'
fmax = 'auto'
colormap = 'gist_rainbow'

# Define plot and font size 
rc('text',usetex=True)
rc('font',size=font)

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

print sim.q
print sim.m_box
x = np.zeros([nx])
y = np.zeros([ny])

#------------------------------------------------------------------------
# Some setup 
#
for i in range(nx):
   x[i] = i*2*np.pi/(nx-1)
for j in range(ny):
   y[j] = j*2*np.pi/(ny-1)

xp = np.zeros([nx,ny])
yp = np.zeros([nx,ny])
zp = np.zeros([nx,ny])

for i in range(nx):
   for j in range(ny):
      xp[i,j] = x[i]/(2*np.pi)*sim.length
      yp[i,j] = y[j]/np.abs(sim.ky[1])
      zp[i,j] = 0.0

aspect = np.max(yp)/np.max(xp)

# Generate vector of time frames 
if istr == '-1':
    ivec = range(nt)
else:
    ivec = str2list(istr)

u=specmap(sim.mass[species],sim.z[species])

# Set filename root and title
isfield = True
if (moment == 'n'):
    fdata = '.cgyro.kxky_n'
    title = r'${\delta \mathrm{n}}_'+u+'$'
    isfield = False
elif (moment == 'e'):
    fdata = '.cgyro.kxky_e'
    title = r'${\delta \mathrm{E}}_'+u+'$'
    isfield = False
elif (moment == 'phi'):
    fdata = '.cgyro.kxky_phi'
    title = r'$\delta\phi$'
elif (moment == 'apar'):
    fdata = '.cgyro.kxky_apar'
    title = r'$\delta A_\parallel$'
elif (moment == 'bpar'):
    fdata = '.cgyro.kxky_bpar'
    title = r'$\delta B_\parallel$'

# ERROR CHECKS

# Check to see if data exists (try binary data first)
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print 'INFO: (plot_viz) Found binary data in '+fdata 
else:
    print 'ERROR: (plot_viz) No data for -moment '+moment+' exists.  Try -moment phi'
    sys.exit()

if isfield:
    n_chunk = 2*nr*nn*nth
else:
    n_chunk = 2*nr*ns*nn*nth
   
# This is the logic to generate a frame
t = 0.0
def frame():
   
   if isfield:
      a = np.reshape(aa,(2,nr,nth,nn),order='F')
   else:
      a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')

   mlab.figure(size=(800,600))
   if isfield:
      c = a[0,:,:,:]+1j*a[1,:,:,:]
   else:
      c = a[0,:,:,species,:]+1j*a[1,:,:,species,:]
                
   f = np.zeros([nx,ny],order='F')
   
   gapy.torcut(c,f)
          
   if fmin == 'auto':
      f0=np.min(f)
      f1=np.max(f)
   else:
      f0=float(fmin)
      f1=float(fmax)

   mlab.mesh(xp,yp,zp,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=0.4)
   print 'INFO: (plot_torcut) min=%e , max=%e  (t=%e)' % (f0,f1,t)

   if ftype == 'screen':
      mlab.show()
   else:
      fname = fdata+str(i)
      # Filename uses frame number 
      plt.savefig(str(i)+'.'+ftype)
      # Close each time to prevent memory accumulation
      plt.close()
                
i = 0

work = True
# Open binary file
fbin = open(fdata,'rb') 

while work:
   try:
      aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
   except:
      sys.exit()
      
   i = i+1
   print 'INFO: (plot_viz) Time index '+str(i) 
   if i in ivec:
      frame()
   if i == max(ivec):
      sys.exit()
