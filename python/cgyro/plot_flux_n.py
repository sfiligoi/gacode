import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w = float(sys.argv[2])
moment = sys.argv[3]

sim = cgyrodata('./')

ns = sim.n_species

fig = plt.figure(figsize=(6*ns,6))

color = ['m','k','b','c']

ky  = sim.ky

dk  = ky[1]-ky[0]
ave = np.zeros((sim.n_n))

if moment == 'n':
    imoment = 0 
    mtag = '\Gamma'
    y = sim.flux_n
elif moment == 'e':
    imoment = 1
    mtag = 'Q'
    y = sim.flux_e
elif moment == 'm':
    print 'm not implemented.'
    sys.exit()
elif moment == 's':
    print 's not implemented.'
    sys.exit()
else:
    print 'ERROR (plot_flux_n.py) Invalid moment.'
    sys.exit()

# Determine tmin
imin=0
for i in range(len(sim.t)):
    if sim.t[i] < (1.0-w)*sim.t[len(sim.t)-1]:
        imin = i+1

for ispec in range(ns):
    stag = str(ispec)
    ax = fig.add_subplot(1,ns,ispec+1)
    ax.set_xlabel(r'$k_\theta \rho_s$')
    ax.set_ylabel(r'$'+mtag+'_'+stag+'$',color='k')
    for j in range(sim.n_n):
        ave[j] = average(y[ispec,j,:],sim.t,w)
    ax.set_title(r'$'+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'$')
    ax.bar(ky-dk/2.0,ave,width=dk/1.1,color=color[ispec],alpha=0.5,edgecolor='black')
        
if ftype == 'screen':
   plt.show()
else:
   outfile = 'flux_n.'+ftype
   plt.savefig(outfile)
