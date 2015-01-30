import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data import GYROData

sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])
window    = float(sys.argv[4])
plotfile  = sys.argv[5]
datafile  = sys.argv[6]

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
ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k',fontsize=GFONTSIZE)
#ax.set_title(r'$\mathtt{'+sys.argv[1]+'}$')
#=====================================

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i

color = ['k','m','b','c']

# Loop over species
if datafile == 'none':
    # Plot data to screen or image file.
    for i in range(n_kinetic):
        ave   = average(flux0[i,i_moment,:],t,window)
        stag  = sim.tagspec[i]
        label = stag+': '+str(round(ave,3))
        y     = ave*np.ones(len(t))
        ax.plot(t[imin:],y[imin:],'--',color=color[i])
        ax.plot(t,flux0[i,i_moment,:],label=label,color=color[i])
else:
    # Write data to datafile
    print 'INFO: (gyro_plot) Output to datafile not supported.  Use raw out.gyro.gbflux.'

ax.legend()
if plotfile == 'screen':
    plt.show()
else:
    plt.savefig(plotfile)
    print "INFO: (gyro_plot) Wrote plot to " + plotfile + "."
