import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
itime = int(sys.argv[2])
ispec = int(sys.argv[3])
tmax = float(sys.argv[4])

sim = cgyrodata('./')

if itime > sim.n_time-1:
    itime = sim.n_time-1

fig = plt.figure(figsize=(12,12))
fig.subplots_adjust(left=0.07,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)
fig.suptitle(r'${\rm species}='+str(ispec)+'$')

func = sim.hb

p = 0
for row in range(3):

    p = p+1

    if row == 0:
        ie = 0
        ix = 0
    if row == 1:
        ie = sim.n_energy/2
        ix = sim.n_xi/2
    if row == 2:
        ie = sim.n_energy-1
        ix = sim.n_xi-1

    #========================================================
    ax = fig.add_subplot(3,3,p)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    ax.set_title(r'$\xi=0 \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\theta/\pi$')

    if sim.n_xi%2 == 0:
        hp = np.array(func[:,:,ispec,sim.n_xi/2,ie,itime]+
                      func[:,:,ispec,sim.n_xi/2-1,ie,itime])*0.5
    else:
        hp = np.array(func[:,:,ispec,sim.n_xi/2,ie,itime])
 
    ax.plot(sim.thetab/np.pi,hp[0,:],'-o',color='black',markersize=2)
    ax.plot(sim.thetab/np.pi,hp[1,:],'-o',color='blue',markersize=2)

    if sim.n_radial > 1:
        if tmax < 0.0:
            ax.set_xlim([1-sim.n_radial,-1+sim.n_radial])
        else:
            ax.set_xlim([-tmax,tmax])
    else:
        ax.set_xlim([1,3])

    #========================================================

    p = p+1

    #========================================================
    ax = fig.add_subplot(3,3,p)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    ax.set_title(r'$\theta=0 \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\xi = v_\parallel/v$')

    n0 = (sim.n_radial/2)*sim.n_theta+sim.n_theta/2 

    hp = np.array(func[0,:,ispec,:,ie,itime])
    ax.plot(sim.xi,hp[n0,:],'-o',color='black',markersize=2)
    hp = np.array(func[1,:,ispec,:,ie,itime])
    ax.plot(sim.xi,hp[n0,:],'-o',color='blue',markersize=2)
    ax.set_xlim([-1,1])
    #========================================================

    p = p+1

    #========================================================
    ax = fig.add_subplot(3,3,p)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    ax.set_title(r'$\theta=0 \quad {\rm ix}='+str(ix)+'$')
    ax.set_xlabel(r'$x=\sqrt{\varepsilon}$')

    n0 = (sim.n_radial/2)*sim.n_theta+sim.n_theta/2

    hp = np.array(func[0,:,ispec,ix,:,itime])
    ax.plot(np.sqrt(sim.energy),hp[n0,:],'-o',color='black',markersize=2)
    hp = np.array(func[1,:,ispec,ix,:,itime])
    ax.plot(np.sqrt(sim.energy),hp[n0,:],'-o',color='blue',markersize=2)
    #========================================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'hbcut.'+ftype
    plt.savefig(outfile)
