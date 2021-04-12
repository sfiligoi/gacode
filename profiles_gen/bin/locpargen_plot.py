import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text',usetex=True)
rc('font',size=16)

fig = plt.figure(figsize=(8,6)) ; ax = fig.add_subplot(111)

data = np.genfromtxt('out.locpargen.theta')
x =data[:,0]/np.pi
rr=data[:,1]
rt=data[:,2]
zr=data[:,3]
zt=data[:,4]

ax.set_xlabel(r'$\theta$')
ax.set_ylabel(r'$F$')
ax.grid(which="both",ls=":")
ax.grid(which="major",ls=":")

ax.plot(x,rr*zt-zr*rt,color='k')
ax.plot(x,rr*zt,color='r')
ax.plot(x,-zr*rt,color='b')
ax.plot(x,x*0.0,linestyle='--')

#ax.legend(loc=2)

ax.set_xlim([-1,1])
#ax.set_ylim([0.1,2.2])

plt.tight_layout()
plt.show()

