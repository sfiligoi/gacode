import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text',usetex=True)
rc('font',size=16)

if len(sys.argv) == 1:
    ftype = 'screen'
else:
    ftype = sys.argv[1]

#=====================================
fig = plt.figure(figsize=(12,6))
fig.suptitle(r'$r/a=0.56\;\mathrm{with\;collisions}$')
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$k_\theta \rho_s$')
ax.set_ylabel(r'$(a/c_s)\,\gamma$')
ax2.grid(which="majorminor",ls=":")
ax2.grid(which="major",ls=":")
ax2.set_xlabel(r'$k_\theta \rho_s$')
ax2.set_ylabel(r'$(a/c_s)\,\omega$')
fig.subplots_adjust(left=0.07,right=0.96,top=0.9,bottom=0.14)
#=====================================

itg=r'$\mathrm{ITG}$'
tem=r'$\mathrm{TEM}$'


root='1itg'
data = np.loadtxt(root+"/fieldeigen_param.out")
k = data[:]
data = np.loadtxt(root+"/fieldeigen_omega.out")
w = data[:,0]
g = data[:,1]

ax.plot(k,g,'-m',label=itg,linewidth=2,alpha=0.4)
ax2.plot(k,w,'-m',linewidth=2,alpha=0.4)

root='1itga'
data = np.loadtxt(root+"/fieldeigen_param.out")
k = data[:]
data = np.loadtxt(root+"/fieldeigen_omega.out")
w = data[:,0]
g = data[:,1]

ax.plot(k,g,'-m',linewidth=2,alpha=0.4)
ax2.plot(k,w,'-m',linewidth=2,alpha=0.4)

root='2tem'
data = np.loadtxt(root+"/fieldeigen_param.out")
k = data[:]
data = np.loadtxt(root+"/fieldeigen_omega.out")
w = data[:,0]
g = data[:,1]

ax.plot(k,g,'-b',label=tem,linewidth=2,alpha=0.4)
ax2.plot(k,w,'-b',linewidth=2,alpha=0.4)

ax.set_xlim([0,4])
ax.set_ylim([0,0.5])
ax2.set_xlim([0,4])
ax2.set_ylim([-0.5,1.1])

ax.legend(loc=2)

if ftype == 'screen':
    plt.show()
else:
    outfile = 'plot.'+ftype
    plt.savefig(outfile)


