import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w = float(sys.argv[2])
moment = sys.argv[3]
ymax = sys.argv[4]

sim = cgyrodata('./')
sim.getbigflux()

ns = sim.n_species

fig = plt.figure(figsize=(6*ns,6))

color = ['m','k','b','c']

ky  = sim.ky
ave = np.zeros((sim.n_n,ns))

if moment == 'n':
    imoment = 0 
    mtag = '\Gamma'
    y = np.sum(sim.kxky_flux_n,axis=0)
    fname = 'out.cgyro.ky_flux_n'
elif moment == 'e':
    imoment = 1
    mtag = 'Q'
    y = np.sum(sim.kxky_flux_e,axis=0)
    fname = 'out.cgyro.ky_flux_e'
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

if ftype == 'dump':
    # Datafile output

    arr = np.zeros([len(ky),ns+1])
    arr[:,0] = ky
    stag = '# (k_y rho_s'
    for ispec in range(ns):
        for j in range(sim.n_n):
            ave = average(y[ispec,j,:],sim.t,w)
            arr[j,ispec+1] = ave
        stag = stag+' , s'+str(ispec)
            
    fid = open(fname,'w')
    fid.write('# Moment  : '+mtag+'\n')
    fid.write('# Time    : '+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'\n')
    fid.write(stag+')\n')
    np.savetxt(fid,arr,fmt='%.5e')
    fid.close()

    print 'INFO: (plot_ky_flux) Wrote output to '+fname
    sys.exit()
else:

    if ky[1] < 0.0:
        ky = -ky
        xlabel=r'$-k_\theta \rho_s$'
    else:
        xlabel=r'$k_\theta \rho_s$'
    
    dk = ky[1]-ky[0]
    
    for ispec in range(ns):
        for j in range(sim.n_n):
            ave[j,ispec] = average(y[ispec,j,:],sim.t,w)

    for ispec in range(ns):
        stag = str(ispec)
        ax = fig.add_subplot(1,ns,ispec+1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$'+mtag+'_'+stag+'$',color='k')
        ax.set_title(r'$'+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'$')
        ax.bar(ky-dk/2.0,ave[:,ispec],width=dk/1.1,color=color[ispec],alpha=0.5,edgecolor='black')
        if ymax != 'auto':
            ax.set_ylim([0,float(ymax)])


if ftype == 'screen':
   plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_ky_flux) Created '+fname
    plt.savefig(fname)
