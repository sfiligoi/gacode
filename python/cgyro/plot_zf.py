import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]
ifield = sys.argv[2]

sim = cgyrodata('./')
sim.getbig()

phic = sim.phi[0,0,0,:]

y     = phic[:]/phic[0]
y_ave = np.average(y)

print 'INFO: (plot_zf) Average =',y_ave

#======================================
fig = plt.figure(figsize=(10,6))
fig.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.12)

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$\Phi/\Phi_0$',fontsize=GFONTSIZE)

ax.plot(sim.t,y,color='k')
ax.plot([0,max(sim.t)],[y_ave,y_ave],color='r',label=r'$\mathrm{Average}$')

y_th = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
ax.plot([0,max(sim.t)],[y_th,y_th],color='b',label=r'$\mathrm{RH \; theory}$')

ax.legend()
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'zf.'+str(ifield)+'.'+ftype
    plt.savefig(outfile)

