import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from prgen_fshape import *

# Function to compute Fourier integrals
# f,w are periodic
def moment(n,f,w,d):

    s0 = 0.0
    s1 = 0.0
    for i in range(n-1):
      s0 = s0+0.5*(f[i]*w[i]+f[i+1]*w[i+1])*d[i]
      s1 = s1+0.5*(w[i]*w[i]+w[i+1]*w[i+1])*d[i]

    return s0/s1

def plot(r,z,x,vr,xr,cr,sr):

    nf = len(cr)-1
    
    rmin = xr[0] 
    rmaj = xr[1]
    zmin = xr[2]*rmin
    zmaj = xr[3] 

    # Angles approximated by Fourier series
    n_arc = len(x)
    pr = np.ones(n_arc)*cr[0]
    for p in range(1,nf+1):
        pr = pr+cr[p]*np.cos(p*x)+sr[p]*np.sin(p*x)

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
    # Parameterized contour
    rp = rmaj+rmin*np.cos(x+pr)
    zp = zmaj+zmin*np.sin(x)
    ax.plot(rp,zp,'-m',linewidth=1)

    # PLOT angle 
    ax = fig.add_subplot(122)
    ax.set_xlabel(r'$x = \theta_Z/(2\pi)$')
    ax.set_ylabel(r'$\theta_R-\theta_Z$')
    ax.grid(which="both",ls=":")

    x=x/(2*np.pi)

    ax.set_xlim([0,1])

    ax.plot(x,pr,'-m',linewidth=1)

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

    plt.tight_layout()
    plt.show()


def oldfourier(ri,zi,nf,rnorm):

   npsi = len(rnorm)
   
   ari = np.zeros([nf+1,npsi])
   bri = np.zeros([nf+1,npsi])
   azi = np.zeros([nf+1,npsi])
   bzi = np.zeros([nf+1,npsi])

   for i in range(npsi-1):
      r=ri[:,i+1] ; z=zi[:,i+1]
      ar,br,az,bz = prgen_fshape(r,z,nf)
      ari[:,i+1] = ar[:]
      bri[:,i+1] = br[:]
      azi[:,i+1] = az[:]
      bzi[:,i+1] = bz[:]

   # Repair origin
   ari[0,:] = extrap(rnorm,ari[0,:]) 
   azi[0,:] = extrap(rnorm,azi[0,:]) 
   for i in range(1,nf+1):
      ari[i,:] = zero(rnorm,ari[i,:]) 
      bri[i,:] = zero(rnorm,bri[i,:]) 
      azi[i,:] = zero(rnorm,azi[i,:]) 
      bzi[i,:] = zero(rnorm,bzi[i,:]) 

   u = ari
   u = np.append(u,bri)
   u = np.append(u,azi)
   u = np.append(u,bzi)
   u.tofile('fluxfit.geo')

   if 0==1:
      fig = plt.figure(figsize=(12,8))

      # PLOT contour
      ax = fig.add_subplot(111)
      for i in range(1,4):
         ax.plot(rnorm,ari[i,:],'-k',linewidth=1)
         ax.plot(rnorm,azi[i,:],'-k',linewidth=1)
         ax.plot(rnorm,bri[i,:],'-k',linewidth=1)
         ax.plot(rnorm,bzi[i,:],'-k',linewidth=1)
      plt.show()
   
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
