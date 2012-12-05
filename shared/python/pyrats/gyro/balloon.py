"""This file is executed by the bash script gyro_plot when a plot of the
fluxes is requested."""

import sys
import numpy as np
from gacodeplotdefs import *
from pyrats.gyro.data import GYROData

sim    = GYROData(sys.argv[1])
index  = int(sys.argv[2])
ftype  = sys.argv[3]

sim.read_balloon()
print sim.balloon.keys()

key = sim.balloon.keys()[index]

ytitle = '\\begin{verbatim}'+key+'\\end{verbatim}'

if key == 'balloon_a':
    ytitle = r'$\delta A_\parallel$'
if key == 'balloon_phi':
    ytitle = r'$\delta \phi$'
if key == 'balloon_aperp':
    ytitle = r'$\delta B_\parallel$'
if key == 'balloon_epar':
    ytitle = r'$\delta E_\parallel$'

#======================================
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$')
ax.set_ylabel(ytitle)
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
