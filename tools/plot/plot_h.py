import sys
import numpy as np
from gacodeplotdefs import *
import matplotlib.cm as cm

ftype = sys.argv[1]
itime = int(sys.argv[2])
ispec = int(sys.argv[3])

#-------------------------------------------------------
# Read grid dimension and axes
#
data = np.loadtxt('out.cgyro.grids')

n_species = int(data[0])
n_radial  = int(data[1])
n_theta   = int(data[2])
n_energy  = int(data[3])
n_xi      = int(data[4])

indx_r = np.array(data[5:5+n_radial],dtype=int)

mark  = 5+n_radial
theta = np.array(data[mark:mark+n_theta])

mark   = mark+n_theta
energy = np.array(data[mark:mark+n_energy])

mark = mark+n_energy
xi   = np.array(data[mark:mark+n_xi])

mark = mark+n_xi
thetab = np.array(data[mark:mark+n_theta*n_radial])
#-------------------------------------------------------


#-------------------------------------------------------
# Read H
#
data = np.loadtxt('out.cgyro.hx')
hx = np.reshape(data,(2,n_radial*n_theta,n_species,n_xi,n_energy),'F')
#-------------------------------------------------------

fig = plt.figure(figsize=(14,12))
fig.subplots_adjust(left=0.1,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)
fig.suptitle(r'${\rm species}='+str(ispec)+'$')

p = 0
for row in range(3):

    p = p+1

    if row == 0:
        ie = 0
    if row == 1:
        ie = n_energy/2
    if row == 2:
        ie = n_energy-1

    #======================================
    ax = fig.add_subplot(3,2,p)

    ax.set_title(r'${\rm Re}(h) \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\theta/\pi$')
    ax.set_ylabel(r'$\xi = v_\parallel/v$')

    hp = np.transpose(np.array(hx[0,:,ispec,:,ie]))
    hmin = hp.min()
    hmax = hp.max()
    dh = (hmax-hmin)/100.0

    levels = np.arange(hmin-dh,hmax+dh,dh)

    ax.contourf(thetab/np.pi,xi,hp,levels,cmap=cm.jet,origin='lower')
    ax.set_xlim([1-n_radial,-1+n_radial])

    # Plot dots for mesh points
    if (row == 1):
        for i in range(n_theta*n_radial):
            for j in range(n_xi):
                ax.plot([thetab[i]/np.pi],[xi[j]],marker='.',color='k',markersize=4)

    #======================================

    p = p+1

    #======================================
    ax = fig.add_subplot(3,2,p)

    
    ax.set_title(r'${\rm Im}(h) \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\theta/\pi$')
    ax.set_ylabel(r'$\xi = v_\parallel/v$')

    hp = np.transpose(np.array(hx[1,:,ispec,:,ie]))
    hmin = hp.min()
    hmax = hp.max()
    dh = (hmax-hmin)/100.0

    levels = np.arange(hmin-dh,hmax+dh,dh)

    ax.contourf(thetab/np.pi,xi,hp,levels,cmap=cm.jet,origin='lower')
    ax.set_xlim([1-n_radial,-1+n_radial])
    #======================================


if ftype == 'screen':
    plt.show()
else:
    outfile = 'h.'+ftype
    plt.savefig(outfile)
