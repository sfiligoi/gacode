import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]

sim = cgyrodata('./')

#======================================
# Set figure size and axes
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.09,right=0.96,top=0.92,bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$')
ax.set_ylabel(r'$\mathrm{Integration~Error}$')
ax.set_yscale('log')
#======================================

ax.plot(sim.t[2:],sim.error[2:])
ax.set_xlim([0,max(sim.t)])

if ftype == 'screen':
    plt.show()
else:
    outfile = 'out.cgyro.error.'+ftype
    plt.savefig(outfile)
