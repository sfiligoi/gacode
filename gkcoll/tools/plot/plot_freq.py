import sys
import numpy as np

import matplotlib.pyplot as plt

import math

GFONTSIZE=18
ftype  = sys.argv[1]

#-------------------------------------------------------
# Read grid dimension and axes
#
data = np.loadtxt('out.gkcoll.freq')
#-------------------------------------------------------

#-------------------------------------------------------
# Read time
t = np.loadtxt('out.gkcoll.time')
n_time = len(t)
#-------------------------------------------------------

fig = plt.figure(figsize=(12,6))

#======================================
ax = fig.add_subplot(121)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$(a/c_s)\, \omega$',fontsize=GFONTSIZE)

ax.plot(t,data[:,0],color='k')

ax.set_xlim([0,max(t)])
#======================================
#======================================
ax = fig.add_subplot(122)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$(a/c_s)\, \gamma$',fontsize=GFONTSIZE)

ax.plot(t,data[:,1],color='k')

ax.set_xlim([0,max(t)])
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = key+'.'+ftype
    plt.savefig(outfile)
