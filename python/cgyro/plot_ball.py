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
    ytag = r'$\delta\phi$'
elif ifield == 1:
    f = sim.aparb[0,:,itime]+1j*sim.aparb[1,:,itime]
    ytag = r'$\delta A_\parallel$'
elif ifield == 2:
    f = sim.bparb[0,:,itime]+1j*sim.bparb[1,:,itime]
    ytag = r'$\delta B_\parallel$'

fig = plt.figure(figsize=(10,6))

#======================================
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.13,right=0.96,top=0.92,bottom=0.12)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$\theta_*/\pi$')
ax.set_ylabel(ytag)


if sim.n_radial == 1:
    # Manage n=0 (ZF) case
    x = sim.theta/np.pi
    ax.set_xlim([-1,1])
else:
    # Assume n > 0 (ballooning mode) if n_radial > 1
    x = sim.thetab/np.pi
    if tmax < 0.0:
        ax.set_xlim([1-sim.n_radial,-1+sim.n_radial])
    else:
        ax.set_xlim([-tmax,tmax])

y1 = np.real(f)
y2 = np.imag(f)

ax.plot(x,y1,'-o',color='black',markersize=2,label=r'$\mathrm{Re}$')
ax.plot(x,y2,'-o',color='red',markersize=2,label=r'$\mathrm{Im}$')

ax.legend()
#======================================

if ftype == 'screen':
    plt.show()
elif ftype == 'dump':
    data = np.column_stack((x,y1,y2))
    np.savetxt('out.cgyro.dump',data,fmt='%.8e')
    print 'INFO: (plot_ball) Created out.cgyro.dump'
else:
    outfile = 'out.cgyro.ball.'+str(ifield)+'.'+ftype
    print 'INFO: (plot_ball) Created '+outfile
    plt.savefig(outfile)
