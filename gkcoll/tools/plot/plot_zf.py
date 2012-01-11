import sys
import numpy as np

import matplotlib.pyplot as plt

import math

GFONTSIZE=18
ftype  = sys.argv[1]

#-------------------------------------------------------
# Read grid dimension and axes
#
data = np.loadtxt('out.gkcoll.grids')

n_species = int(data[0])
n_radial  = int(data[1])
n_theta   = int(data[2])
n_energy  = int(data[3])
n_xi      = int(data[4])

indx_r = np.array(data[5:5+n_radial],dtype=int)

mark  = 5+n_radial
theta = np.array(data[mark:mark+n_theta])

mark   = mark+n_theta
energy = np.array(data[mark:mark+n_energy])

mark = mark+n_energy
xi   = np.array(data[mark:mark+n_xi])

mark = mark+n_xi
thetab = np.array(data[mark:mark+n_theta*n_radial])
#-------------------------------------------------------

#-------------------------------------------------------
# Read time
t = np.loadtxt('out.gkcoll.time')
n_time = len(t)

#-------------------------------------------------------

#-------------------------------------------------------
# Read phib
#
data = np.loadtxt('out.gkcoll.phi')
phib = np.reshape(data,(2,n_theta,n_radial,n_time),'F')
# Construct complex eigenfunction at selected time
phic = phib[0,n_theta/4,0,:]

#-------------------------------------------------------

fig = plt.figure(figsize=(6,6))

#======================================
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$',fontsize=GFONTSIZE)

ax.plot(t,np.real(phic),color='k',label='Re')

#ax.set_xlim([1-n_radial,-1+n_radial])
ax.legend()
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = key+'.'+ftype
    plt.savefig(outfile)
