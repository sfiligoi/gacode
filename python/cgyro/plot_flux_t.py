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

color = ['k','m','b','c','g','r']

if moment == 'n':
    mtag = '\Gamma'
    ttag = 'G'
    y = np.sum(sim.flux_n,axis=(0,2))
elif moment == 'e':
    mtag = 'Q'
    ttag = 'Q'
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

ax.set_title(r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$')

for ispec in range(ns):
    ave = average(y[ispec,:],t,w)
    y_ave = ave*np.ones(len(t))
    label = r'$'+mtag+'_'+str(ispec)+'/'+mtag+'_\mathrm{GB}: '+str(round(ave,3))+'$'
    # Average
    ax.plot(t[imin:],y_ave[imin:],'--',color=color[ispec])
    # Time trace
    ax.plot(sim.t,y[ispec,:],label=label,color=color[ispec])

ax.legend(loc=2)

if ymax != 'auto':
    ax.set_ylim([0,float(ymax)])

fname = 'out.cgyro.fluxt.'+ttag+'.'

if ftype == 'screen':
   plt.show()
elif ftype == 'dump':
    data = np.column_stack((sim.t,y[0,:]))
    head = '(cs/a) t     '+ttag+'_1/'+ttag+'_GB'
    for ispec in range(1,ns,1):
        head = head+'       '+ttag+'_'+str(ispec+1)+'/'+ttag+'_GB'
        data = np.column_stack((data,y[ispec,:]))
    np.savetxt('out.cgyro.dump.flux_t',data,fmt='%.8e',header=head)
    print 'INFO: (plot_flux_t) Created '+fname+'txt'
else:
   outfile = fname+ftype
   print 'INFO: (plot_flux_t) Created '+outfile
   plt.savefig(outfile)
