import sys
import numpy as np
import matplotlib.cm as cm

import matplotlib.pyplot as plt

GFONTSIZE=18
ftype = sys.argv[1]
itime = int(sys.argv[2])
ispec = int(sys.argv[3])


#-------------------------------------------------------
# Read grid dimension and axes
#
data = np.loadtxt('out.gkcoll.grids')

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
data = np.loadtxt('out.gkcoll.hx')
hx = np.reshape(data,(2,n_radial*n_theta,n_xi,n_energy,n_species),'F')
#-------------------------------------------------------

fig = plt.figure(figsize=(12,12),dpi=100)

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
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    if row == 0:
        ax.set_title(r'$\xi=0 \quad {\rm spec}='+str(ispec)+'$')
    ax.set_xlabel(r'$\theta/\pi$')

    hp = np.array(hx[0,:,:,ie,ispec])
    ax.plot(thetab,hp[:,n_xi/2],'-o',color='black',markersize=4)
    hp = np.array(hx[1,:,:,ie,ispec])
    ax.plot(thetab,hp[:,n_xi/2],'-o',color='blue',markersize=4)
    ax.set_xlim([1-n_radial,-1+n_radial])
    #======================================

    p = p+1

    #======================================
    ax = fig.add_subplot(3,2,p)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    if row == 0:
         ax.set_title(r'$\theta=0 \quad {\rm spec}='+str(ispec)+'$')
    ax.set_xlabel(r'$\xi = v_\parallel/v$')

    n0 = (n_radial/2)*n_theta+n_theta/2

    hp = np.array(hx[0,:,:,ie,ispec])
    ax.plot(xi,hp[n0,:],'-o',color='black',markersize=4)
    hp = np.array(hx[1,:,:,ie,ispec])
    ax.plot(xi,hp[n0,:],'-o',color='blue',markersize=4)
    ax.set_xlim([-1,1])
    #======================================

# Cut 1

if ftype == 'screen':
    plt.show()
else:
    outfile = 'h.'+ftype
    plt.savefig(outfile)
