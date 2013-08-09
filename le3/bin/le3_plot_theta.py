import sys
import os
import numpy as np
from gacodeplotdefs import *

simdir  = sys.argv[1]
imgfile = sys.argv[2]

rc('lines',linewidth=1)
#rc('mathtext',default='bf')

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
print mt,mp

for j in range(mp):
    ax.plot(t[:]/np.pi,(tb[:,j]-t[:])/np.pi)

TICKS=[0,1,2]
LABELS=[r'\boldmath{$0$}',r'\boldmath{$\pi$}',r'\boldmath{$2\pi$}']
ax.set_xticks(TICKS)
ax.set_xticklabels(LABELS)
ax.set_xlim(0,2)

TICKS=[-0.125,0,0.125]
LABELS=[r'\boldmath{$\pi/8$}',r'\boldmath{$0$}',r'\boldmath{$\pi/8$}']
ax.set_yticks(TICKS)
ax.set_yticklabels(LABELS)
ax.set_ylim(-0.125,0.125)

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print "INFO: (le3_plot_theta) Wrote plot to "+imgfile+"."
