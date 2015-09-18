import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
field = sys.argv[2]

sim = cgyrodata('./')
sim.getbig()

fig = plt.figure(figsize=(10,8))

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$k_y$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$k_x$',fontsize=GFONTSIZE)

# Note array structure
# self.phi = np.reshape(data,(2,self.n_radial,self.n_n,nt),'F')

nx=sim.n_radial
ny=sim.n_n

#f = sim.phisq[:,:,-1]
f = sim.flux_e[:,0,:,-1]

#f[nx/2,0] = 1e-12
#f[0,0]    = 1e-12
#f = np.log10(f)
fmax = f.max()
fmin = 0.0

d = (fmax-fmin)/100.0
levels = np.arange(fmin-d,fmax+d,d)

ax.contourf(sim.ky,sim.kx,f,levels,origin='lower')

if ftype == 'screen':
    plt.show()
else:
    outfile = 'kxky_phisq.'+ftype
    plt.savefig(outfile)
