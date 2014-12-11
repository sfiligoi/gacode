import sys
import numpy as np
from gacodeplotdefs import *

ftype = sys.argv[1]

#-------------------------------------------------------
# Read grid dimension and axes
#
data = np.loadtxt('out.cgyro.grids')

n_species = int(data[0])
n_radial  = int(data[1])
n_theta   = int(data[2])
n_energy  = int(data[3])
n_xi      = int(data[4])

indx_r = np.array(data[5:5+n_radial],dtype=int)

mark  = 5+n_radial
theta = np.array(data[mark:mark+n_theta])

mark   = mark+n_theta
energy = np.array(data[mark:mark+n_energy])

mark = mark+n_energy
xi   = np.array(data[mark:mark+n_xi])

mark = mark+n_xi
thetab = np.array(data[mark:mark+n_theta*n_radial])
#-------------------------------------------------------

#-------------------------------------------------------
# Read time
t = np.loadtxt('out.cgyro.time')
n_time = len(t)

#-------------------------------------------------------

#-------------------------------------------------------
# Read phib
#
data = np.loadtxt('out.cgyro.phi')
phib = np.reshape(data,(2,n_theta,n_radial,n_time),'F')
# Construct complex eigenfunction at selected time
phic = phib[0,n_theta/4,0,:]
y = np.real(phic)/np.real(phic[0])
y_ave = np.sum(y)/len(y)
print y_ave

y_th = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
#-------------------------------------------------------

fig = plt.figure(figsize=(10,6))

#======================================
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)

ax.plot(t,y,color='k',label='Re')
ax.plot([0,max(t)],[y_ave,y_ave],color='r',label='ave')
#ax.plot([0,max(t)],[y_th,y_th],color='b',label='RH theory')

ax.legend()
#======================================

#for i in range(len(y)):
#    print y[i]

if ftype == 'screen':
    plt.show()
else:
    outfile = key+'.'+ftype
    plt.savefig(outfile)

