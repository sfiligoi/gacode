import sys
import os
import numpy as np
import matplotlib.cm as cm
from gacodeplotdefs import *

simdir  = sys.argv[1]
imgfile = sys.argv[2]
index   = int(sys.argv[3])

rc('lines',linewidth=1)

dim = np.loadtxt(simdir+'/out.le3.geoscalar')
vec = np.loadtxt(simdir+'/out.le3.geovector')

nts = int(dim[0])
nps = int(dim[1])
ns  = int(dim[2])
q   = float(dim[7])

if index == 0:
    symbol = '\\bar{\\theta}-\\theta'
if index == 1:
    symbol = 'vx'
if index == 2:
    symbol = 'flux'
if index == 3:
    symbol = 'B'
if index == 4:
    symbol = '1'
if index == 5:
    symbol = 'Ave'
if index == 6:
    symbol = 'B'

z = np.zeros(ns)
z = vec[:,index]

nx = 256

t = 2*np.pi*np.arange(nx)/float(nx-1)

dx = 2*np.ones(nts+1)
dy = 2*np.ones(nps+1)
dx[0] = 1
dy[0] = 1

# Plotting
fig = plt.figure(figsize=(9,6))

ax = fig.add_subplot(111)

for s in range(4):
    p = q*(t+2*np.pi*s)
    f = np.zeros(nx)

    i = 0
    for ips in range(nps+1):
        for its in range(nts+1):
            if its > 0:
                # A
                f = f+np.sin(its*t)*np.cos(ips*p)*z[i]*dx[its]*dy[ips]
                i = i+1
            if ips > 0 and its > 0:
                # B
                f = f+np.sin(its*t)*np.sin(ips*p)*z[i]*dx[its]*dy[ips]
                i = i+1
            if ips+its > 0:
                # C
                f = f+np.cos(its*t)*np.cos(ips*p)*z[i]*dx[its]*dy[ips]
                i = i+1
            if ips > 0:
                # D
                f = f+np.cos(its*t)*np.sin(ips*p)*z[i]*dx[its]*dy[ips]
                i = i+1
            if ips+its == 0:
                f = f+z[i]
                i = i+1


    ax.plot(t/np.pi,f)

ax.set_xlabel(r'$\theta/\pi$')
ax.set_xlabel(r'$'+symbol+'$')
ax.set_xlim([0,2])

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print "INFO: (le3_plot_fieldline) Wrote plot to "+imgfile+"."
