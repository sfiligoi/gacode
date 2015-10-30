import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data import GYROData

sim   = GYROData(sys.argv[1])
w     = float(sys.argv[2])
ftype = sys.argv[3]

# Read freq data
sim.read_moment_u()

ntheta = sim.profile['n_theta_plot']
nx = sim.profile['n_x']

print nx,ntheta
y = np.real(sim.moment_u[ntheta/3,nx/3,0,0,:])

y = y[:]/y[0]
t = sim.t['(c_s/a)t']

#----------------------------------------------------
# Average calculations

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-w)*t[len(t)-1]:
        imin = i

ave = average(y[:],t,w)
print 'INFO: (plot_zf) Integral time-average = %.6f' % ave

ave_vec = ave*np.ones(len(t))
#----------------------------------------------------

#======================================
fig = plt.figure(figsize=(10,6))
fig.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.12)

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$\Phi/\Phi_0$',fontsize=GFONTSIZE)

ax.plot(t,y,color='k')
ax.plot(t[imin:],ave_vec[imin:],color='b',label=r'$\mathrm{Average}$',linewidth=1)

theory = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
ax.plot([0,max(t)],[theory,theory],color='grey',label=r'$\mathrm{RH \; theory}$',alpha=0.3,linewidth=4)

ax.legend()
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'zf.'+ftype
    plt.savefig(outfile)

