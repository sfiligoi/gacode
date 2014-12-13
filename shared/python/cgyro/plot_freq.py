import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]

sim = cgyrodata('./')

fig = plt.figure(figsize=(12,6))

#======================================
# Omega
ax = fig.add_subplot(121)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$(a/c_s)\, \omega$',fontsize=GFONTSIZE)

ax.plot(sim.t,sim.freq[:,0],color='k')

ax.set_xlim([0,max(sim.t)])
#======================================

#======================================
# Gamma
ax = fig.add_subplot(122)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$(a/c_s)\, \gamma$',fontsize=GFONTSIZE)

ax.plot(sim.t,sim.freq[:,1],color='k')

ax.set_xlim([0,max(sim.t)])
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'freq.'+ftype
    plt.savefig(outfile)
