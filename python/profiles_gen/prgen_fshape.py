import sys
import string
import numpy as np

def prgen_fshape(rd,zd,nf):

    nd = len(rd)

    s = np.argmax(rd)

    # Shift elements so that first index at max(R).
    rd[0:-1] = np.roll(rd[0:-1],-s) ; rd[-1] = rd[0] 
    zd[0:-1] = np.roll(zd[0:-1],-s) ; zd[-1] = zd[0]
    
    # Construct equally-spaced poloidal angle
    theta  = np.linspace(0,1,nd)*2*np.pi
    dtheta = theta[1]-theta[0]
    
    ar = np.zeros(nf+1)
    br = np.zeros(nf+1)
    az = np.zeros(nf+1)
    bz = np.zeros(nf+1)

    ds = dtheta/np.pi

    # Trapezoidal integration spectrally-accurate
    for n in range(nf+1):
       y = np.cos(n*theta[:])*rd[:] ; ar[n] = ds*np.sum(y[:-1])
       y = np.sin(n*theta[:])*rd[:] ; br[n] = ds*np.sum(y[:-1])
       y = np.cos(n*theta[:])*zd[:] ; az[n] = ds*np.sum(y[:-1])
       y = np.sin(n*theta[:])*zd[:] ; bz[n] = ds*np.sum(y[:-1])

    return ar,br,az,bz

