import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
field = sys.argv[2]

sim = cgyrodata('./')

fig = plt.figure(figsize=(12,6))

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$ky$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$p$',fontsize=GFONTSIZE)

p = np.arange(sim.n_radial)

f = sim.pwr_phi[:,:,-1]

fmin = f.min()
fmax = f.max()
d = (fmax-fmin)/100.0
levels = np.arange(fmin-d,fmax+d,d)

ax.contourf(sim.ky,p,f,levels,cmap=cm.jet,origin='lower')

if ftype == 'screen':
    plt.show()
else:
    outfile = 'pwr_np.'+ftype
    plt.savefig(outfile)
