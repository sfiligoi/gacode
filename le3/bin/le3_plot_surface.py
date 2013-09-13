import sys
import os
import numpy as np
import matplotlib.cm as cm
from gacodeplotdefs import *

simdir  = sys.argv[1]
imgfile = sys.argv[2]

xm = np.loadtxt(simdir+'/out.le3.t')
ym = np.loadtxt(simdir+'/out.le3.p')
zm = np.loadtxt(simdir+'/out.le3.b')

x = xm/np.pi
y = ym/np.pi

nx = len(x)
ny = len(y)

z = zm.reshape((nx,ny),order='F')

fig = plt.figure(figsize=(6,6))
fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)

ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.set_title(r'$B(\theta,\varphi)$',size=18)

dz = 0.005
levels = np.arange(np.amin(z)-dz,np.amax(z)+dz,dz)
ax.contourf(x,y,z,levels,cmap=cm.jet,origin='lower')
ax.set_xlabel(r'$\varphi/\pi$')
ax.set_ylabel(r'$\theta/\pi$')
ax.set_xlim([0,2])
ax.set_ylim([0,2])

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig('b3d.eps')
    print "INFO: (le3_plot_surface) Wrote plot to "+imgfile+"."

#-------------------------------------------------------------------

zm = np.loadtxt(simdir+'/out.le3.tb')

z = zm.reshape((nx,ny),order='F')

fig = plt.figure(figsize=(6,6))
fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)

ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.set_title(r'$\bar{\theta}(\theta,\varphi)$',size=18)

dz = 0.005
levels = np.arange(np.amin(z)-dz,np.amax(z)+dz,dz)
ax.contourf(x,y,z,levels,cmap=cm.jet,origin='lower')
ax.set_xlabel(r'$\varphi/\pi$')
ax.set_ylabel(r'$\theta/\pi$')
ax.set_xlim([0,2])
ax.set_ylim([0,2])

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig('t3d.eps')
    print "INFO: (le3_plot_surface) Wrote plot to "+imgfile+"."
