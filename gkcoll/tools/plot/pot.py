"""This file is executed by the bash script gyro_plot when a plot of the
fluxes is requested."""

import sys
import numpy as np

import matplotlib.pyplot as plt

GFONTSIZE=18
ftype = 'screen'

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
data = np.loadtxt('out.gkcoll.phiB')
phib = np.reshape(data,(2,n_theta*n_radial,n_time),'F')
#-------------------------------------------------------

#-------------------------------------------------------
# Read H
#
data = np.loadtxt('out.gkcoll.hx')
hx = np.reshape(data,(2,n_species,n_radial,n_theta,n_energy,n_xi),'F')
#-------------------------------------------------------


fig = plt.figure(figsize=(12,12))

#======================================
ax = fig.add_subplot(221)
ax.set_title('p=0')
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$',fontsize=GFONTSIZE)

ax.plot(thetab/np.pi,phib[0,:,-1],color='k',label='Re')
ax.plot(thetab/np.pi,phib[1,:,-1],color='c',label='Im')

ax.set_xlim([1-n_radial,-1+n_radial])
ax.legend()
#======================================
#======================================
ax = fig.add_subplot(222)
ax.set_title('Real, p=-1')
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$',fontsize=GFONTSIZE)

for i in range(n_xi):
    ax.plot(theta/np.pi,hx[0,0,2,:,4,i])

ax.set_xlim([-1,1])
#ax.legend()
#======================================
#======================================
ax = fig.add_subplot(223)
ax.set_title('Im, p=0')
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$',fontsize=GFONTSIZE)

for i in range(n_xi):
    ax.plot(theta/np.pi,hx[0,0,3,:,4,i])

ax.set_xlim([-1,1])
#ax.legend()
#======================================
#======================================
ax = fig.add_subplot(224)
ax.set_title('Re, p=1')
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$',fontsize=GFONTSIZE)

for i in range(n_xi):
    ax.plot(theta/np.pi,hx[0,0,4,:,4,i])

ax.set_xlim([-1,1])
ax.legend()
#======================================

#if ftype == 'screen':
plt.show()
#else:
#    outfile = key+'.'+ftype
#    plt.savefig(outfile)
