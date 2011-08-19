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

GFONTSIZE=16

sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])
window    = float(sys.argv[4])

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])
n_n       = int(sim.profile['n_n'])

t    = sim.t['(cbar_s/a)t']
flux = sim.gbflux_n

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag = '\mathrm{total}'
else:
    i_field = int(field)
    flux0 = flux[:,i_field,:,:,:]
    if i_field == 0: 
        ftag = '\mathrm{electrostatic}'
    if i_field == 1: 
        ftag = '\mathrm{flutter}'
    if i_field == 2: 
        ftag = '\mathrm{compression}'

# Manage moment
if i_moment == 0: 
    mtag = '\Gamma/\Gamma_\mathrm{GB}'
if i_moment == 1: 
    mtag = 'Q/Q_\mathrm{GB}'
if i_moment == 2: 
    mtag = '\Pi/\Pi_\mathrm{GB}'
if i_moment == 3: 
    mtag = 'S/S_\mathrm{GB}'

#======================================
fig = plt.figure(figsize=(6*n_kinetic,6))
#=====================================

color = ['k','m','b','c']
k = sim.profile['kt_rho']
dk = k[1]-k[0]

# Loop over species
for i in range(n_kinetic):
    ax = fig.add_subplot(1,n_kinetic,i+1)
    ax.set_xlabel(r'$k_\theta \rho_s$',fontsize=GFONTSIZE)
    ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k',fontsize=GFONTSIZE)
    if i == n_kinetic-1:
        stag = 'elec  '
    else:
        stag = 'ion-'+str(i+1)+' '

    ave = np.zeros((n_n))
    for j in range(n_n):
        ave[j] = average(flux0[i,i_moment,j,:],t,window)

    label = stag
    ax.set_title(stag)
    ax.bar(k-dk/2.0,ave,width=dk/1.1,color=color[i],alpha=0.4,edgecolor='black')

#ax.legend()
plt.show()
