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

sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])
window    = float(sys.argv[4])
ftype     = sys.argv[5]

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])

t = sim.t['(c_s/a)t']

# Read data in gbflux_i and make gbflux
sim.read_gbflux_i()
sim.make_gbflux()

flux = sim.gbflux

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag  = sim.tagfield[3]
else:
    i_field = int(field)
    flux0   = flux[:,i_field,:,:]
    ftag    = sim.tagfield[i_field]

# Manage moment
mtag = sim.tagmom[i_moment]

#======================================
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel('$(c_s/a) t$',fontsize=GFONTSIZE)
ax.set_ylabel('$'+mtag+' \;('+ftag+')$',color='k',fontsize=GFONTSIZE)
ax.set_title(sys.argv[1])
#=====================================

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i

color = ['k','m','b','c']

# Loop over species
for i in range(n_kinetic):
    ave   = average(flux0[i,i_moment,:],t,window)
    stag  = sim.tagspec[i]
    label = stag+': '+str(round(ave,3))
    y     = ave*np.ones(len(t))
    ax.plot(t[imin:],y[imin:],'--',color=color[i])
    ax.plot(t,flux0[i,i_moment,:],label=label,color=color[i])

ax.legend()
if ftype == 'screen':
    plt.show()
else:
    outfile = 'gbflux.'+ftype
    plt.savefig(outfile)
