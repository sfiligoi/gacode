import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]

sim = cgyrodata('./')
sim.getgeo()

fig = plt.figure(figsize=(15,9))
fig.subplots_adjust(left=0.07,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)

theta = sim.geo[:,0]/np.pi

for p in range(8):
    p1 = p+1
    y = sim.geo[:,p1]

    ax = fig.add_subplot(2,4,p1)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")
    ax.set_xlabel(r'$\theta/\pi$',fontsize=GFONTSIZE)
    ax.set_title(r'$'+sim.geotag[p1]+'$')
    ax.plot(theta,y)
    ax.set_xlim([-1,1])

#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'geo.'+ftype
    plt.savefig(outfile)
