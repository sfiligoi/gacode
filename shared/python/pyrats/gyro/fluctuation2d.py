"""This file is executed by the bash script gyro_plot when a plot of the
fluxes is requested."""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.gyro.data import GYROData
import scipy.interpolate

GFONTSIZE=18

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

if field == 's':
    field = 0

field = int(field)

if n_theta_plot == 1: 
    j = 0
else:
    j = n_theta_plot/2

# Read data in moment_u
if i_moment > 2:
    # Plot fields
    sim.read_moment_u()
    # f[n_x,n_field,n_n]
    f = np.array(sim.moment_u[j,:,:,:,-1],dtype=complex)
else:
    # f[n_x,n_kinetic,n_n]
    if i_moment == 0:
        sim.read_moment_n()
        f = np.array(sim.moment_n[j,:,:,:,-1],dtype=complex)
    if i_moment == 1:
        sim.read_moment_v()
        f = np.array(sim.moment_v[j,:,:,:,-1],dtype=complex)
    if i_moment == 2:
        sim.read_moment_e()
        f = np.array(sim.moment_e[j,:,:,:,-1],dtype=complex)

# Plot resolution
nxp = 512
nyp = 16*n_n

# Plot (x,y) arrays
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

#======================================
if i_moment > 2:
    n_plot = n_field
else:
    n_plot = n_kinetic

fig = plt.figure(figsize=(7*n_plot,6))

for ip in range(n_plot):

    for i_n in np.arange(n_n):
        fr = scipy.interpolate.splrep(sim.profile['r'],np.real(f[:,ip,i_n]))
        fi = scipy.interpolate.splrep(sim.profile['r'],np.imag(f[:,ip,i_n]))
        ffine[i_n,:] = scipy.interpolate.splev(xp,fr)+1j*scipy.interpolate.splev(xp,fi)

        for i in np.arange(nxp):
            for j in np.arange(nyp):
                zp[i,j] = np.sum(np.real(ffine[:,i]*phase[:,j]))

    ax = fig.add_subplot(1,n_plot,ip+1)
    ax.pcolorfast(xp,yp,np.transpose(zp),cmap=plt.cm.jet)
    ax.set_xlim([np.min(xp),np.max(xp)])
    ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)
    ax.set_ylabel(r'$\alpha/2\pi$',fontsize=GFONTSIZE)
   
#=====================================



#ax.legend()
if ftype == 'screen':
    plt.show()
else:
    outfile = 'fluctuation2d.'+ftype
    plt.savefig(outfile)
