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

