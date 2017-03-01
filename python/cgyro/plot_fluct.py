import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
moment = sys.argv[2]
ymin = sys.argv[3]
ymax = sys.argv[4]
nx = int(sys.argv[5])
ny = int(sys.argv[6])

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species

#(2,n_radial,n_species,n_n,nt)
n_chunk = 2*nr*ns*nn
n = 0
i = 0

i_time = 65
i_spec = 0

a = np.zeros([n_chunk])
for line in open('out.cgyro.kxky_n'):
    a[n] = float(line)
    n = n+1
    if n == n_chunk:
        i = i+1
        n = 0
        if i == i_time:
            break

a = np.reshape(a,(2,nr,ns,nn),'F')
c = a[0,:,i_spec,:]+1j*a[1,:,i_spec,:]
#c[nr,nn]

f = np.zeros([nx,ny])
epx = np.zeros([nx,nr],dtype=np.complex64)
eny = np.zeros([ny,nn],dtype=np.complex64)

for i in range(nx):
    x = i*2*np.pi/(nx-1)
    for p in range(nr):    
        epx[i,p]=np.exp(1j*(p-nr/2)*x)

for j in range(ny):
    y = j*2*np.pi/(ny-1)
    for n in range(nn):    
        eny[j,n]=np.exp(-1j*n*y)

for j in range(ny):
    for p in range(nr):
        for n in range(nn):
            f[:,j] = f[:,j]+np.real(c[p,n]*epx[:,p]*eny[j,n])
 
fig = plt.figure(figsize=(15,15))
fig.subplots_adjust(left=0.05,right=0.96,top=0.91,bottom=0.14)
ax = fig.add_subplot(111)
ax.set_xlabel(r'$x/2\pi$',fontsize=16)
ax.set_ylabel(r'$y/2\pi$',fontsize=16)

f0=np.min(f)
f1=np.max(f)

levels = np.arange(f0,f1,(f1-f0)/128)
ax.contourf(np.transpose(f),levels,cmap=plt.cm.jet)

plt.show()
