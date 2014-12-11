import sys
import numpy as np
from gacodeplotdefs import *

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
m_box     = int(data[5])

indx_r = np.array(data[6:6+n_radial],dtype=int)

mark  = 6+n_radial
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
hx = hx/np.max(hx)
#-------------------------------------------------------

fig = plt.figure(figsize=(12,12))
fig.subplots_adjust(left=0.07,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)
fig.suptitle(r'${\rm species}='+str(ispec)+'$')

p = 0
for row in range(3):

    p = p+1

    if row == 0:
        ie = 0
        ix = 0
    if row == 1:
        ie = n_energy/2
        ix = n_xi/2
    if row == 2:
        ie = n_energy-1
        ix = n_xi-1

    #======================================
    ax = fig.add_subplot(3,3,p)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    ax.set_title(r'$\xi=0 \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\theta/\pi$')

    if n_xi%2 == 0:
        hp = np.array(hx[:,:,ispec,n_xi/2,ie]+hx[:,:,ispec,n_xi/2-1,ie])*0.5
    else:
        hp = np.array(hx[:,:,ispec,n_xi/2,ie])
 
    ax.plot(thetab/np.pi,hp[0,:],'-o',color='black',markersize=2)
    ax.plot(thetab/np.pi,hp[1,:],'-o',color='blue',markersize=2)

    if n_radial > 1:
        ax.set_xlim([1-n_radial,-1+n_radial])
    else:
        ax.set_xlim([1,3])

    #======================================

    p = p+1

    #======================================
    ax = fig.add_subplot(3,3,p)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    ax.set_title(r'$\theta=0 \quad {\rm ie}='+str(ie)+'$')
    ax.set_xlabel(r'$\xi = v_\parallel/v$')

    n0 = (n_radial/2)*n_theta+n_theta/2 

    hp = np.array(hx[0,:,ispec,:,ie])
    ax.plot(xi,hp[n0,:],'-o',color='black',markersize=2)
    hp = np.array(hx[1,:,ispec,:,ie])
    ax.plot(xi,hp[n0,:],'-o',color='blue',markersize=2)
    ax.set_xlim([-1,1])
    #======================================

    p = p+1

    #======================================
    ax = fig.add_subplot(3,3,p)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")

    ax.set_title(r'$\theta=0 \quad {\rm ix}='+str(ix)+'$')
    ax.set_xlabel(r'$x=\sqrt{\varepsilon}$')

    n0 = (n_radial/2)*n_theta+n_theta/2

    hp = np.array(hx[0,:,ispec,ix,:])
    ax.plot(np.sqrt(energy),hp[n0,:],'-o',color='black',markersize=2)
    #ax.plot(np.sqrt(energy),hp[n0,:]*np.exp(-energy),'-o',color='black',markersize=2)
    hp = np.array(hx[1,:,ispec,ix,:])
    ax.plot(np.sqrt(energy),hp[n0,:],'-o',color='blue',markersize=2)
    #ax.plot(np.sqrt(energy),hp[n0,:]*np.exp(-energy),'-o',color='blue',markersize=2)
    #======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'h.'+ftype
    plt.savefig(outfile)
