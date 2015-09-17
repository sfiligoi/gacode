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
ax.set_xlabel(r'$ky$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$p$',fontsize=GFONTSIZE)

# Note array structure
# self.phi = np.reshape(data,(2,self.n_theta,self.n_radial,self.n_n,nt),'F')

f = sim.pwr_phi[:,:,-1]
fmax = f.max()

# Zero wavenumebrs
f[0,0] = 1.0
f[sim.n_radial/2,0] = 1.0

f = np.log(f)

fmin = f.min()
fmax = np.log(fmax)

d = (fmax-fmin)/100.0
levels = np.arange(fmin-d,fmax+d,d)

ax.contourf(sim.ky,sim.kx,f,levels,origin='lower')

if ftype == 'screen':
    plt.show()
else:
    outfile = 'pwr_np.'+ftype
    plt.savefig(outfile)
