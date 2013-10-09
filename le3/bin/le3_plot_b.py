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

# vec[:,0] = thetabar-theta
# vec[:,1] = vdriftx
# vec[:,2] = flux
# vec[:,3] = uparB
# vec[:,4] = upar 
# vec[:,5] = fsa
# vec[:,6] = bmag

z = np.zeros(ns)
z = vec[:,1]

nplot = 64

t = 2*np.pi*np.arange(nplot)/float(nplot-1)
p = 2*np.pi*np.arange(nplot)/float(nplot-1)

i = 0
f = np.zeros([nplot,nplot])

for ips in range(nps+1):
    for its in range(nts+1):
        if its > 0:
            # A
            f = f+2*np.sin(its*t)*np.cos(ips*p)*z[i]
            i = i+1
        if ips > 0 and its > 0:
            # B
            f = f+2*np.sin(its*t)*np.sin(ips*p)*z[i]
            i = i+1
        if ips+its > 0:
            # C
            f = f+2*np.cos(its*t)*np.cos(ips*p)*z[i]
            i = i+1
        if ips > 0:
            # D
            f = f+2*np.cos(its*t)*np.sin(ips*p)*z[i]
            i = i+1
        if ips+its == 0:
            f = f+z[i]
            i = i+1


# phi=0
for i in range(nplot):
    print t[i]/np.pi,f[0,i]

# Plotting
fig = plt.figure(figsize=(6,6))
fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)

ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.set_title(r'$B(\theta,\varphi)$',size=18)

dz=(np.amax(f)-np.amin(f))/1000
levels = np.arange(np.amin(f)-dz,np.amax(f)+dz,dz)
ax.contourf(t/np.pi,p/np.pi,f,levels,cmap=cm.jet,origin='lower')
ax.set_ylabel(r'$\varphi/\pi$')
ax.set_xlabel(r'$\theta/\pi$')
ax.set_xlim([0,2])
ax.set_ylim([0,2])

plt.show()
