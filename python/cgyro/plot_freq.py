import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]

sim = cgyrodata('./')

fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.085,right=0.97,top=0.92,bottom=0.12)

#======================================
# Omega
ax = fig.add_subplot(121)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(TIME)
ax.set_ylabel(r'$(a/c_s)\, \omega$')

for i in range(sim.n_n):
    ax.plot(sim.t,sim.freq[0,i,:])

ax.set_xlim([0,max(sim.t)])
#======================================

#======================================
# Gamma
ax = fig.add_subplot(122)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(TIME)
ax.set_ylabel(r'$(a/c_s)\, \gamma$')

for i in range(sim.n_n):
    ax.plot(sim.t,sim.freq[1,i,:])

ax.set_xlim([0,max(sim.t)])
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'out.cgyro.freq.'+ftype
    plt.savefig(outfile)
