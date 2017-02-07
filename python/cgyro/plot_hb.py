import sys
import numpy as np
from gacodeplotdefs import *
import matplotlib.cm as cm
from cgyro.data import cgyrodata

ftype = sys.argv[1]
itime = int(sys.argv[2])
ispec = int(sys.argv[3])
tmax = float(sys.argv[4])
mesh = int(sys.argv[5])

theta=0.0

sim = cgyrodata('./')

if itime > sim.n_time-1:
    itime = sim.n_time-1

fig = plt.figure(figsize=(14,12))
fig.subplots_adjust(left=0.1,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)
fig.suptitle(r'${\rm species}='+str(ispec)+'$')

# Compute index for theta value in pitch angle and energy plots
i0 = int(round((1.0+float(theta))*sim.n_theta/2.0))
if i0 > sim.n_theta-1:
    i0 = sim.n_theta-1
n0 = (sim.n_radial/2)*sim.n_theta+i0

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

    hp = np.transpose(np.array(sim.hb[0,:,ispec,:,ie,itime]))
    h_norm = 0.5*(hp[sim.n_xi/2-1,n0]+hp[sim.n_xi/2,n0])
    hp = hp/h_norm
    hmin = hp.min()
    hmax = hp.max()
    dh = (hmax-hmin)/100.0

    levels = np.arange(hmin-dh,hmax+dh,dh)

    ax.contourf(sim.thetab/np.pi,sim.xi,hp,levels,cmap=cm.jet,origin='lower')
    if tmax < 0.0:
        ax.set_xlim([1-sim.n_radial,-1+sim.n_radial])
    else:
        ax.set_xlim([-tmax,tmax])

    # Plot dots for mesh points
    if row == 1 and mesh == 1:
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

    hp = np.transpose(np.array(sim.hb[1,:,ispec,:,ie,itime]))
    hmin = hp.min()
    hmax = hp.max()
    dh = (hmax-hmin)/100.0

    levels = np.arange(hmin-dh,hmax+dh,dh)

    ax.contourf(sim.thetab/np.pi,sim.xi,hp,levels,cmap=cm.jet,origin='lower')
    if tmax < 0.0:
        ax.set_xlim([1-sim.n_radial,-1+sim.n_radial])
    else:
        ax.set_xlim([-tmax,tmax])
    #======================================


if ftype == 'screen':
    plt.show()
else:
    outfile = 'hb.'+ftype
    plt.savefig(outfile)
