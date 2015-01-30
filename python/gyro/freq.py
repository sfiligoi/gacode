import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data import GYROData

fig = plt.figure(figsize=(16,8))

GFONTSIZE=16

sim       = GYROData(sys.argv[1])
window    = float(sys.argv[2])
ftype     = sys.argv[3]

t    = sim.t['(c_s/a)t']

# Read freq data
sim.read_freq()


# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i

color = ['k','m','b','c']

#======================================
ax = fig.add_subplot(121)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$(a/c_s)\gamma$',color='k',fontsize=GFONTSIZE)
#=====================================

#Gamma
for i in range(sim.profile['n_n']):
    ax.plot(t[imin:],sim.freq['(a/c_s)gamma'][i,imin:],color=color[i])

#======================================
ax = fig.add_subplot(122)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$(a/c_s)\omega$',color='k',fontsize=GFONTSIZE)
#=====================================

#Omega
for i in range(sim.profile['n_n']):
    ax.plot(t[imin:],sim.freq['(a/c_s)w'][i,imin:],color=color[i])

if ftype == 'screen':
    plt.show()
else:
    outfile = 'freq.'+ftype
    plt.savefig(outfile)
