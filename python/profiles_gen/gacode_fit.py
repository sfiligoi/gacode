import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Function to compute Fourier integrals
# f,w are periodic 
def moment(n,f,w,d):
    
    s0 = 0.0
    s1 = 0.0
    for i in range(n-1):
      s0 = s0+0.5*(f[i]*w[i]+f[i+1]*w[i+1])*d[i]
      s1 = s1+0.5*(w[i]*w[i]+w[i+1]*w[i+1])*d[i]

    return s0/s1

def fit(r,z,n_arc,nf,pflag):
   
   # Number of theta-points for plotting
   dz = np.zeros(n_arc)
   ur = np.zeros(n_arc) ; uz = np.zeros(n_arc)
   vr = np.zeros(n_arc) ; vz = np.zeros(n_arc)

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
      for i in range(n_arc):
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
   xr[2] = zmin
   xr[3] = zmaj

   for i in range(n_arc-1):
      dz[i] = vz[i+1]-vz[i]

   dz[-1] = dz[0]

   # Compute expansion coefficients (cr):
   #  vr = sum cr*cos(nx)+sr*sin(nx)
   cr = np.zeros(nf+1) 
   sr = np.zeros(nf+1) 

   cr[0] = moment(n_arc,vr,np.ones(n_arc),dz)
   for p in range(1,nf+1):
      cr[p] = moment(n_arc,vr,np.cos(p*x),dz)
      sr[p] = moment(n_arc,vr,np.sin(p*x),dz)

   # Regenerate based on shape param
   pr = np.ones(n_arc)*cr[0]
   for p in range(1,nf+1):
      pr = pr+cr[p]*np.cos(p*x)+sr[p]*np.sin(p*x)

   # Parameterized contour
   rp = rmaj+rmin*np.cos(x+pr)
   zp = zmaj+zmin*np.sin(x)

   if pflag:
      #--------------------------------------------------------------------
      # PLOTTING
      #--------------------------------------------------------------------
      # Latex fonts
      rc('text',usetex=True)
      rc('font',size=18)

      fig = plt.figure(figsize=(18,9))

      # PLOT contour
      ax = fig.add_subplot(121,aspect='equal')
      ax.set_xlabel(r'$R$')
      ax.set_ylabel(r'$Z$')
      ax.grid(which="both",ls=":")

      # Data
      ax.plot(r,z,'-k',linewidth=2,alpha=0.3)
      # Parameterized
      ax.plot(rp,zp,'-m',linewidth=1)

      # PLOT angle 
      ax = fig.add_subplot(122)
      ax.set_xlabel(r'$x = \theta_Z/(2\pi)$')
      ax.set_ylabel(r'$\theta_R-\theta_Z$')
      ax.grid(which="both",ls=":")

      x=x/(2*np.pi)

      ax.set_xlim([0,1])

      # Shift elements so that x<1.0
      s=-2
      for i in range(n_arc):
         if x[i] > 1.0:
            s = i-1
            break

      if s > -2:
         x[0:-1]  = np.roll(x[0:-1],-s)  ; x[-1] = x[0] 
         vr[0:-1] = np.roll(vr[0:-1],-s) ; vr[-1] = vr[0]
         pr[0:-1] = np.roll(pr[0:-1],-s) ; pr[-1] = pr[0]

      for i in range(n_arc-1):
         if x[n_arc-2-i] > x[n_arc-1-i]:
            x[n_arc-2-i] = x[n_arc-2-i]-1.0

      # Angles from data
      ax.plot(x,vr,'-k',linewidth=2,alpha=0.3)
      # Angles approximated by Fourier series
      ax.plot(x,pr,'-m',linewidth=1)

      plt.tight_layout()
      plt.show()

   return cr,sr,xr
