#----------------------------------------------------------------------
# gacodeplotdefs.py
#
# PURPOSE:
#  Simplify/standardize gacode plotting.
#----------------------------------------------------------------------

import matplotlib.pyplot as plt
from matplotlib import rc

GFONTSIZE=18
rc('text',usetex=True)
rc('font',size=GFONTSIZE)

#---------------------------------------------------------------
def average(f,t,window):
 
    n_time = len(t)

    # Manage case with 2 time points (eigenvalue)
    if len(t) == 2:
        tmin = t[n_time-1]
        tmax = tmin
        ave  = f[n_time-1]
        return ave

    tmin = (1.0-window)*t[n_time-1]
    tmax = t[n_time-1]

    t_window = 0.0
    ave      = 0.0
    for i in range(n_time-1):
        if t[i] > tmin: 
            ave = ave+0.5*(f[i]+f[i+1])*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window

    return ave
#---------------------------------------------------------------
#---------------------------------------------------------------
def average_n(f,t,window,n):
 
    import numpy as np

    ave = np.zeros(n)

    n_time = len(t)
    tmin = (1.0-window)*t[n_time-1]
    tmax = t[n_time-1]

    t_window = 0.0
    for i in range(n_time-1):
        if t[i] > tmin:
            ave[:] = ave[:]+0.5*(f[:,i]+f[:,i+1])*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window

    return ave
#---------------------------------------------------------------
#---------------------------------------------------------------
def smooth_pro(x,z,p,n):

    import numpy as np

    nx = len(x)
    xf = np.zeros((nx-1)*n+1)
    zf = np.zeros((nx-1)*n+1)
    pf = np.zeros((nx-1)*n+1)
    j = 0
    for i in range(nx-1):
        for m in range(n):
            u = m/(1.0*n)
            xf[j] = x[i]*(1-u)+x[i+1]*u
            zf[j] = z[i]*(1-u)+z[i+1]*u
            j = j+1

    xf[j] = x[nx-1]
    zf[j] = z[nx-1]
    pf[j] = p[nx-1]

    # Exponential integration to obtain smooth profiles
    for i in np.arange(j,0,-1):
        pf[i-1] = pf[i]*np.exp(0.5*(xf[i]-xf[i-1])*(zf[i]+zf[i-1]))


    return xf,pf
#---------------------------------------------------------------
