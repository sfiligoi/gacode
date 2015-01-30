import sys
import numpy as np
from gacodeplotdefs import *
import matplotlib.cm as cm
from cgyro.data import cgyrodata

ftype = sys.argv[1]
itime = int(sys.argv[2])
ispec = int(sys.argv[3])

sim = cgyrodata('./')

fig = plt.figure(figsize=(14,12))
fig.subplots_adjust(left=0.1,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)
fig.suptitle(r'${\rm species}='+str(ispec)+'$')

p = 0
for row in range(3):

    p = p+1

    if row == 0:
        ie = 0
    if row == 1:
        ie = sim.n_energy/2
    if row == 2:
        ie = sim.n_energy-1

    #======================================
    ax = fig.add_subplot(3,2,p)

    ax.set_title(r'${\rm Re}(h) \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\theta/\pi$')
    ax.set_ylabel(r'$\xi = v_\parallel/v$')

    hp = np.transpose(np.array(sim.hb[0,:,ispec,:,ie]))
    hmin = hp.min()
    hmax = hp.max()
    dh = (hmax-hmin)/100.0

    levels = np.arange(hmin-dh,hmax+dh,dh)

    ax.contourf(sim.thetab/np.pi,sim.xi,hp,levels,cmap=cm.jet,origin='lower')
    ax.set_xlim([1-sim.n_radial,-1+sim.n_radial])

    # Plot dots for mesh points
    if (row == 1):
        for i in range(sim.n_theta*sim.n_radial):
            for j in range(sim.n_xi):
                ax.plot([sim.thetab[i]/np.pi],[sim.xi[j]],marker='.',color='k',markersize=4)

    #======================================

    p = p+1

    #======================================
    ax = fig.add_subplot(3,2,p)

    
    ax.set_title(r'${\rm Im}(h) \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\theta/\pi$')
    ax.set_ylabel(r'$\xi = v_\parallel/v$')

    hp = np.transpose(np.array(sim.hb[1,:,ispec,:,ie]))
    hmin = hp.min()
    hmax = hp.max()
    dh = (hmax-hmin)/100.0

    levels = np.arange(hmin-dh,hmax+dh,dh)

    ax.contourf(sim.thetab/np.pi,sim.xi,hp,levels,cmap=cm.jet,origin='lower')
    ax.set_xlim([1-sim.n_radial,-1+sim.n_radial])
    #======================================


if ftype == 'screen':
    plt.show()
else:
    outfile = 'hb.'+ftype
    plt.savefig(outfile)
