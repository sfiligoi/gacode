import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
field = sys.argv[2]

sim = cgyrodata('./')
sim.getbig()

fig = plt.figure(figsize=(12,6))

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$')
ax.set_ylabel(r'$\left| \Phi \right|$')
ax.set_yscale('log')

y = np.sum(sim.phisq,axis=0)

for n in range(sim.n_n):
    num = r'$n='+str(n)+'$'
    if n==0:
        ax.plot(sim.t,np.sqrt(y[n,:]),linewidth=2,label=num)
    else:
        ax.plot(sim.t,np.sqrt(y[n,:]),label=num)

ax.set_xlim([0,max(sim.t)])
#======================================

if sim.n_n > 16:
    ax.legend(loc=4, ncol=5, prop={'size':12})
else:
    ax.legend(loc=4, ncol=6, prop={'size':12})


if ftype == 'screen':
    plt.show()
else:
    outfile = 'pwr.'+ftype
    plt.savefig(outfile)
