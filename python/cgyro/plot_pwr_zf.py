import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]

sim = cgyrodata('./')
sim.getbig()

fig = plt.figure(figsize=(12,6))

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$')
ax.set_ylabel(r'$\Phi$')
ax.set_yscale('log')

y0 = sim.phisq[:,0,:]

for i in range(sim.n_radial):
    num = r'$p='+str(i)+'$'
    if abs(i-sim.n_radial/2) < 3:
        ax.plot(sim.t,y0[i,:])
#    if i == 0:
#        ax.plot(sim.t,y0[i,:],'r',alpha=0.2)
#    if i == sim.n_radial/2:
#        ax.plot(sim.t,y0[i,:],'--k')

ax.set_xlim([0,max(sim.t)])
#ax.set_ylim([-0.1,0.1])
#======================================

#if sim.n_n > 16:
#    ax.legend(loc=4, ncol=5, prop={'size':12})
#else:
#    ax.legend(loc=4, ncol=6, prop={'size':12})


if ftype == 'screen':
    plt.show()
else:
    outfile = 'pwr_zf.'+ftype
    plt.savefig(outfile)
