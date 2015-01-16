import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]
ifield = sys.argv[2]

sim = cgyrodata('./')
sim.getbig()

phic = sim.phi[0,sim.n_theta/3,0,:]

y     = np.real(phic)/np.real(phic[0])
y_ave = np.sum(y)/len(y)

print len(y)

print 'INFO: (plot_zf) Average =',y_ave

#======================================
fig = plt.figure(figsize=(10,6))

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)

ax.plot(sim.t,y,color='k',label='Re')
ax.plot([0,max(sim.t)],[y_ave,y_ave],color='r',label='ave')

#y_th = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
#ax.plot([0,max(t)],[y_th,y_th],color='b',label='RH theory')

ax.legend()
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'zf.'+str(ifield)+'.'+ftype
    plt.savefig(outfile)

