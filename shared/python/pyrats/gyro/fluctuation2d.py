"""This file is executed by the bash script gyro_plot when a plot of the
fluxes is requested."""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.gyro.data import GYROData
import scipy.interpolate


#---------------------------------------------------------------
def average(f,t,window,n):
 
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

GFONTSIZE=16

sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])
window    = float(sys.argv[4])
ftype     = sys.argv[5]

n_x          = sim.profile['n_x']
n_n          = sim.profile['n_n']
n_field      = sim.profile['n_field']
n_kinetic    = sim.profile['n_kinetic']
n_theta_plot = sim.profile['n_theta_plot']

t = sim.t['(c_s/a)t']

# Read data in gbflux_i
sim.read_moment_u()

if n_theta_plot == 1: 
    j = 0
else:
    j = n_theta_plot/2

f = np.array(sim.moment_u[j,:,0,:,-1],dtype=complex)

nxp = 256
nyp = 8*n_n

xp = np.arange(nxp)/(nxp-1.0)
xp = sim.profile['r'][0]+xp*(sim.profile['r'][-1]-sim.profile['r'][0])

yp    = np.arange(nyp)/(nyp-1.0)
phase = np.zeros((n_n,nyp),dtype=complex)

zp = np.zeros((nxp,nyp))

cn    = 2.0*np.ones(n_n)
cn[0] = 1.0

for i_n in np.arange(n_n):
    phase[i_n,:] = np.exp(-2*np.pi*1j*i_n*yp[:])*cn[i_n]

ffine = np.zeros((n_n,nxp),dtype=complex)

for i_n in np.arange(n_n):
    fr = scipy.interpolate.splrep(sim.profile['r'],np.real(f[:,i_n]))
    fi = scipy.interpolate.splrep(sim.profile['r'],np.imag(f[:,i_n]))
    ffine[i_n,:] = scipy.interpolate.splev(xp,fr)+1j*scipy.interpolate.splev(xp,fi)

for i in np.arange(nxp):
    for j in np.arange(nyp):
        zp[i,j] = np.sum(np.real(ffine[:,i]*phase[:,j]))

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)

#zmin=np.min(zp)
#zmax=np.max(zp)
#dz=(zmax-zmin)/1024
#levels = np.arange(zmin,zmax,dy)
#ax.contourf(re,alpha,np.transpose(y),levels,cmap=plt.cm.jet,origin='lower')
ax.pcolorfast(xp,yp,np.transpose(zp))

ax.set_xlim([np.min(xp),np.max(xp)])
ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$\alpha/2\pi$',fontsize=GFONTSIZE)

# Manage field
if field == 's':
    ftag = '\mathrm{total}'
else:
    i_field = int(field)
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
#ax.grid(which="majorminor",ls=":")
#ax.grid(which="major",ls=":")
#ax.set_title(str(t[imin])+' < (c_s/a) t < '+str(t[-1]))
#=====================================

color = ['k','m','b','c']

# Loop over species
for i in range(n_kinetic):
    if sim.profile['electron_method'] == 2 or  sim.profile['electron_method'] == 4:
        if i == n_kinetic-1:
            stag = 'elec  '
        else:
            stag = 'ion-'+str(i+1)+' '
    if sim.profile['electron_method'] == 1:
        stag = 'ion-'+str(i+1)+' '
    if sim.profile['electron_method'] == 3:
        stag = 'elec  '

    label = stag

#ax.legend()
if ftype == 'screen':
    plt.show()
else:
    outfile = 'fluctuation2d.'+ftype
    plt.savefig(outfile)
