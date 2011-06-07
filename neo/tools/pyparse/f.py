#!/bin/env python

import pylab, numpy as np, mpmath
import matplotlib, matplotlib.cm as cm, matplotlib.mlab as mlab
import matplotlib.pyplot as plt


# Load the data.  Note that pylab loads it by row.  
gridData=pylab.loadtxt("grid.out")
tempfData=pylab.loadtxt("f.out")
# use grid.out info to reshape data
nvars=1
norder=2
nspecies=int(gridData[0])
nenergy=int(gridData[1])
npitch=int(gridData[2]+1)
ntheta=int(gridData[3])
nradial=int(gridData[4+ntheta])
print norder,nradial,nspecies,nenergy,npitch,ntheta
fData=tempfData.reshape(norder,nradial,nspecies,nenergy,npitch,ntheta)
feData=np.zeros((nenergy,npitch), np.float)
fiData=np.zeros((nenergy,npitch), np.float)
nr = 20
nt = 10
for i in range(0,nenergy):
  for j in range(0,npitch):
      fiData[i,j]=fData[0,nr,0,i,j,nt]
      feData[i,j]=fData[0,nr,1,i,j,nt]

radius=np.zeros(nradial, np.float)
for i in range(0,nradial):
  radius[i]=gridData[5+ntheta+i]

theta=np.zeros(ntheta, np.float)
for i in range(0,ntheta):
  theta[i]=gridData[4+i]

# grid from -v_max to v_max in parallel velocity
# and 0 to v_max in perpendicular velocity
energy_max = 16.
v_max = mpmath.sqrt(2.0*energy_max)
delta = v_max/5.
vpll = np.arange(-v_max, v_max+delta/2., delta)
nvpll = len(vpll)
vprp = np.arange(0.0, v_max+delta/2., delta)
nvprp = len(vprp)
fe = np.zeros((nvprp,nvpll), np.float)
fi = np.zeros((nvprp,nvpll), np.float)
f_max = 0.

for n in range(0,nvprp):
  for m in range(0,nvpll):
    v = mpmath.sqrt(vpll[m]**2+vprp[n]**2) 
    if v > v_max:
      fi[n,m] = 0.
      fe[n,m] = 0.
    else:
      xi = vpll[m]/v_max
      u = 2.0*(v/v_max)-1.0
      for i in range(0,nenergy):
        for j in range(0,npitch):
          v_poly = (mpmath.chebyu(i,u)/mpmath.chebyu(i,1.0))*mpmath.legendre(j,xi)
          fm = mpmath.exp(-v**2)
          fi[n,m] = fi[n,m]+fm*fiData[i,j]*v_poly
          fe[n,m] = fe[n,m]+fm*feData[i,j]*v_poly
    print n,m,vpll[m],vprp[n],fi[n,m],fe[n,m]
    
plt.subplot(1,2,1)
CS = plt.contourf(vpll, vprp, fe, 40)
plt.colorbar(CS)
plt.title('Fe')

plt.subplot(1,2,2)
CS = plt.contourf(vpll, vprp, fi, 40)
plt.colorbar(CS)
plt.title('Fi')

plt.show()
