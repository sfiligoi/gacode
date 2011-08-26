"""This file is executed by the bash script gyro_plot when a plot of the
fluxes is requested."""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.gyro.data import GYROData

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

GFONTSIZE=18

sim    = GYROData(sys.argv[1])
index  = int(sys.argv[2])
ftype  = sys.argv[3]

sim.read_balloon()
print sim.balloon.keys()

key = sim.balloon.keys()[index]

#======================================
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$',fontsize=GFONTSIZE)
ax.set_ylabel(key,color='k',fontsize=GFONTSIZE)
#=====================================

n_p   = sim.profile['n_x']/sim.profile['box_multiplier']
n_ang = sim.profile['n_theta_plot']*n_p

x = -(1.0+n_p)+2.0*n_p*np.arange(n_ang)/float(n_ang)

ax.plot(x,np.real(sim.balloon[key][:,0,-1]),color='k',label='Re')
ax.plot(x,np.imag(sim.balloon[key][:,0,-1]),color='m',label='Im')

ax.set_xlim([1-n_p,-1+n_p])

ax.legend()
if ftype == 'screen':
    plt.show()
else:
    outfile = key+'.'+ftype
    plt.savefig(outfile)
