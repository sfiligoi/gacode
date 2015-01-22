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
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$\Phi$',fontsize=GFONTSIZE)

y = np.sum(sim.pwr_phi,axis=0)

for n in range(sim.n_n):
    num = r'$n='+str(n)+'$'
    if n==0:
        ax.plot(sim.t,np.sqrt(y[n,:]),linewidth=2,label=num)
    else:
        ax.plot(sim.t,np.sqrt(y[n,:]),label=num)

ax.set_xlim([0,max(sim.t)])
#======================================

ax.legend(loc=2, ncol=3, fontsize=11, borderpad=2,frameon=False)

if ftype == 'screen':
    plt.show()
else:
    outfile = 'pwr.'+ftype
    plt.savefig(outfile)
