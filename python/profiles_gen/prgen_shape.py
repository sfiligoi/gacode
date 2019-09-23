import sys
import string
import numpy as np

from prgen_shape_util import *

def prgen_shape(r,z,narc,nf,pflag):

   # Number of theta-points for plotting
   dz = np.zeros(narc)
   ur = np.zeros(narc) ; uz = np.zeros(narc)
   vr = np.zeros(narc) ; vz = np.zeros(narc)

   # Pointwise Extrema ; definitions of rmin, rmaj, etc.
   n1 = np.argmax(z) ; m1 = np.argmin(z)
   zmaj = 0.5*(z[n1]+z[m1])
   zmin = 0.5*(z[n1]-z[m1])

   n2 = np.argmax(r) ; m2 = np.argmin(r)
   rmaj = 0.5*(r[n2]+r[m2])
   rmin = 0.5*(r[n2]-r[m2])

   # Index of rightmost point (want this to be 0th index)
   s = n2

   # Shift elements so that first index at max(R).
   r[0:-1] = np.roll(r[0:-1],-s) ; z[0:-1] = np.roll(z[0:-1],-s)
   r[-1] = r[0] ; z[-1] = z[0]

   if z[1] < z[0]:
      # Reverse order (may be needed)
      r = np.flip(r,0) ; z = np.flip(z,0)

   # Compute generalized angles
   eps = 1.0-1e-6
   for i in range(narc):
      # (ur,uz): principle angles (discontinuous)
      uz[i] = np.arcsin(eps*(z[i]-zmaj)/zmin)
      ur[i] = np.arccos(eps*(r[i]-rmaj)/rmin)
      # (vr,vz): proper (continuous) branches
      if i > 4:
         if uz[i] > uz[i-1] and ur[i] > ur[i-1]:
            vz[i] = uz[i] ; vr[i] = ur[i]
         elif uz[i] < uz[i-1] and ur[i] > ur[i-1]:
            vz[i] = np.pi-uz[i] ; vr[i] = ur[i]
         elif uz[i] < uz[i-1] and ur[i] < ur[i-1]:
            vz[i] = np.pi-uz[i] ; vr[i] = 2*np.pi-ur[i]
         elif uz[i] > uz[i-1] and ur[i] < ur[i-1]:
            vz[i] = 2*np.pi+uz[i] ; vr[i] = 2*np.pi-ur[i]
      else:
         vz[i] = uz[i] ; vr[i] = ur[i]

   # Define vr as deviation from vz
   vr = vr-vz ; vr[-1] = vr[0]

   x = vz
   # Now:
   #
   # R = rmaj + rmin*cos(x+vr)
   # Z = zmaj + zmin*sin(x)
   #
   # NOTE: kappa=zmin/rmin

   xr = np.zeros(4)
   xr[0] = rmin
   xr[1] = rmaj
   xr[2] = zmin/rmin
   xr[3] = zmaj

   for i in range(narc-1):
      dz[i] = vz[i+1]-vz[i]

   dz[-1] = dz[0]

   # Compute expansion coefficients (cr):
   #  vr = sum cr*cos(nx)+sr*sin(nx)
   cr = np.zeros(nf+1) 
   sr = np.zeros(nf+1) 

   cr[0] = moment(narc,vr,np.ones(narc),dz)
   for p in range(1,nf+1):
      cr[p] = moment(narc,vr,np.cos(p*x),dz)
      sr[p] = moment(narc,vr,np.sin(p*x),dz)

   if pflag:
       plot(r,z,x,vr,xr,cr,sr)

   return cr,sr,xr
