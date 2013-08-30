import sys
import os
import numpy as np
from gacodeplotdefs import *

simdir  = sys.argv[1]
imgfile = sys.argv[2]

t  = np.loadtxt(simdir+'/out.le3.t')
p  = np.loadtxt(simdir+'/out.le3.p')
tb = np.loadtxt(simdir+'/out.le3.tb')

#======================================
fig = plt.figure(figsize=(7,6))
fig.subplots_adjust(left=0.17, right=0.97, top=0.95, bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'\boldmath{$\theta$}')
ax.set_ylabel(r'\boldmath{$\bar\theta}-\theta$}',color='k')
#=====================================

mt = len(t)
mp = len(p)

x = np.zeros([mt+1])
y = np.zeros([mt+1])

x[0:mt] = t[0:mt]/np.pi
x[mt] = 2.0

for j in range(mp):
    y[0:mt] = (tb[:,j]-t[:])/np.pi
    y[mt] = y[0]
    ax.plot(x,y)

TICKS=[0,1,2]
LABELS=[r'$0$',r'$\pi$',r'$2\pi$']
ax.set_xticks(TICKS)
ax.set_xticklabels(LABELS)
ax.set_xlim(0,2)

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print "INFO: (le3_plot_theta) Wrote plot to "+imgfile+"."
