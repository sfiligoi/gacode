import os,re
import numpy as np

# DMD utility: downsample with averaging (eps is flatness)
def dmddownsample(f,n,eps=0.01):
    npar = f.shape[0]
    nave = f.shape[1]

    # array end trim
    m = (nave//n)*n

    # downsample

    u = np.linspace(-1,1,n)

    w = 1/(1+eps*u*u)
    w = w/np.sum(w)

    fdown = np.sum(f[:,:m].reshape(npar,-1,n)*w[None,None,:],axis=2)

    return fdown
