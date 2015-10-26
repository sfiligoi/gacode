import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w = float(sys.argv[2])
moment = sys.argv[3]
ymax = sys.argv[4]

sim = cgyrodata('./')

ns = sim.n_species
t = sim.t

#======================================
# Set figure size and axes
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$')
#======================================

color = ['m','k','b','c']

if moment == 'n':
    imoment = 0 
    mtag = '\Gamma'
    y = np.sum(sim.flux_n,axis=(0,2))
elif moment == 'e':
    imoment = 1
    mtag = 'Q'
    y = np.sum(sim.flux_e,axis=(0,2))
elif moment == 'm':
    print 'm not implemented.'
    sys.exit()
elif moment == 's':
    print 's not implemented.'
    sys.exit()
else:
    print 'ERROR (plot_flux_t.py) Invalid moment.'
    sys.exit()

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-w)*t[len(t)-1]:
        imin = i

ax.set_title(r'$'+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'$')

for ispec in range(ns):
    ave = average(y[ispec,:],t,w)
    y_ave = ave*np.ones(len(t))
    label = r'$'+mtag+'_'+str(ispec)+': '+str(round(ave,3))+'$'
    # Average
    ax.plot(t[imin:],y_ave[imin:],'--',color=color[ispec])
    # Time trace
    ax.plot(sim.t,y[ispec,:],label=label,color=color[ispec])

ax.legend(loc=2)

if ymax != 'auto':
    ax.set_ylim([0,float(ymax)])

if ftype == 'screen':
   plt.show()
else:
   outfile = 'flux_t.'+ftype
   plt.savefig(outfile)
