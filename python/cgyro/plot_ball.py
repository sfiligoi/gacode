import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]
itime  = int(sys.argv[2])
ifield = int(sys.argv[3])
tmax   = float(sys.argv[4])

sim = cgyrodata('./')

if itime > sim.n_time-1:
    itime = sim.n_time-1

# Construct complex eigenfunction at selected time
if ifield == 0:
    f = sim.phib[0,:,itime]+1j*sim.phib[1,:,itime]
elif ifield == 1:
    f = sim.aparb[0,:,itime]+1j*sim.aparb[1,:,itime]
elif ifield == 2:
    f = sim.bparb[0,:,itime]+1j*sim.bparb[1,:,itime]

fig = plt.figure(figsize=(10,6))

#======================================
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$')

ax.plot(sim.thetab/np.pi,np.real(f),'-o',color='black',markersize=2,label='Re')
ax.plot(sim.thetab/np.pi,np.imag(f),'-o',color='blue',markersize=2,label='Im')

if tmax < 0.0:
    ax.set_xlim([1-sim.n_radial,-1+sim.n_radial])
else:
    ax.set_xlim([-tmax,tmax])

ax.legend()
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'ball.'+str(ifield)+'.'+ftype
    plt.savefig(outfile)
