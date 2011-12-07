"""This file is executed by the bash script gyro_plot when a plot of the
fluxes is requested."""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.gyro.data import GYROData

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

flux = sim.gbflux_i

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag  = sim.tagfield[3]
else:
    i_field = int(field)
    flux0 = flux[:,i_field,:,:,:]
    ftag  = sim.tagfield[i_field]

# Manage moment
mtag = sim.tagmom[i_moment]

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i

#====================================================================
fig = plt.figure(figsize=(7*n_kinetic,6))
#===================================================================

color = ['k','m','b','c']

# Loop over species
for i in range(n_kinetic):
    stag  = sim.tagspec[i]
    ax = fig.add_subplot(1,n_kinetic,i+1)
    ax.set_xlabel(r'$(c_s/a)t$',fontsize=GFONTSIZE)
    ax.set_ylabel(r'$r/a$',fontsize=GFONTSIZE)
    ax.set_title(r'$'+mtag+' \;('+ftag+'\,\mathrm{'+stag+'})$',color='k',fontsize=GFONTSIZE)
    #label = stag+': '+str(round(ave,3))
    #ax.pcolorfast(t[imin:],sim.profile['r'],flux0[i,i_moment,:,imin:],cmap=plt.cm.jet)
    ax.contourf(t[imin:],sim.profile['r'],flux0[i,i_moment,:,imin:],cmap=plt.cm.jet,nlevels=128)
    ax.set_xlim([t[imin],t[-1]])
    ax.set_ylim([sim.profile['r'][0],sim.profile['r'][-1]])

if ftype == 'screen':
    plt.show()
else:
    outfile = 'gbflux.'+ftype
    plt.savefig(outfile)
