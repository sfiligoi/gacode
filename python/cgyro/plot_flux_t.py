import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w = float(sys.argv[2])
moment = sys.argv[3]

sim = cgyrodata('./')

ns = sim.n_species

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$')

color = ['m','k','b','c']

if moment == 'n':
    imoment = 0 
    mtag = '\Gamma'
    y = np.sum(sim.flux_n,axis=1)
elif moment == 'e':
    imoment = 1
    mtag = 'Q'
    y = np.sum(sim.flux_e,axis=1)
elif moment == 'm':
    print 'm not implemented.'
    sys.exit()
elif moment == 's':
    print 's not implemented.'
    sys.exit()
else:
    print 'ERROR (plot_flux_t.py) Invalid moment.'
    sys.exit()

for ispec in range(ns):
    stag = str(ispec)
    ax.plot(sim.t,y[ispec,:],label=r'$'+mtag+'_'+stag+'$')

ax.legend(loc=2, fontsize=12, borderpad=2, frameon=False)

if ftype == 'screen':
   plt.show()
else:
   outfile = 'flux_t.'+ftype
   plt.savefig(outfile)
