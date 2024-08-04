import numpy as np
import sys
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text',usetex=True)
rc('font',size=15)

fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.set_xlabel(r'$T_e$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.subplots_adjust(bottom=0.12,
                    top=0.92,
                    left=0.08,
                    right=0.98)

data = np.loadtxt('he4.txt')
ax.plot(data[:,0],data[:,1],'r')
data = np.loadtxt('He4.txt')
ax.plot(data[:,0],data[:,1],'k')


plt.show()
#plt.savefig(fn+'.png')
#plt.close()

