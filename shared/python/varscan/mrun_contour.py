import sys
import os
import numpy as np
import matplotlib.cm as cm
from gacodeplotdefs import *

simdir  = sys.argv[1]
index   = int(sys.argv[2])
imgfile = sys.argv[3]
code    = sys.argv[4]

if len(sys.argv) > 5:
    x = sys.argv[5]
    var1 = x.split(',')[0]
    var2 = x.split(',')[1] 
else:
    var1 = "x"
    var2 = "y"


nz = 128

xm = np.loadtxt(simdir+'/out.'+code+'.mrun_x')
ym = np.loadtxt(simdir+'/out.'+code+'.mrun_y')
zm = np.loadtxt(simdir+'/out.'+code+'.mrun_z')

x = xm
y = ym

nx = len(x)
ny = len(y)

z = zm[:,index].reshape((ny,nx))

zmin = np.amin(z)
zmax = np.amax(z)
 
dz = (zmax-zmin)/nz
levels = np.arange(zmin,zmax+dz,dz)

fig = plt.figure(figsize=(6,6))
#fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)

ax = fig.add_subplot(111)

f = ax.contourf(x,y,z,levels,cmap=cm.jet,origin='lower')
plt.colorbar(f)
ax.set_xlabel(r'$\mathrm{'+var1+'}$')
ax.set_ylabel(r'$\mathrm{'+var2+'}$')
ax.set_xlim([np.amin(x),np.amax(x)])
ax.set_ylim([np.amin(y),np.amax(y)])

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print "INFO: (mrun_contour) Wrote plot to "+imgfile+"."

#-------------------------------------------------------------------

