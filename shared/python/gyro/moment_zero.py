import sys
import numpy as np
import matplotlib.pyplot as plt
from gyro.data import GYROData

#---------------------------------------------------------------
def average(f,t,window,n):
 
    ave = np.zeros(n)

    n_time = len(t)
    tmin = (1.0-window)*t[n_time-1]
    tmax = t[n_time-1]

    t_window = 0.0
    for i in range(n_time-1):
        if t[i] > tmin: 
            ave[:] = ave[:]+0.5*(f[:,i]+f[:,i+1])*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window

    return ave
#---------------------------------------------------------------

GFONTSIZE=18

sim      = GYROData(sys.argv[1])
i_moment = int(sys.argv[2])
window   = float(sys.argv[3])
ftype    = sys.argv[4]

n_field   = sim.profile['n_field']
n_kinetic = sim.profile['n_kinetic']
n_x       = sim.profile['n_x']

t = sim.t['(c_s/a)t']

# Read data in gbflux_i
sim.read_moment_zero()

delta_n = np.empty([n_x,n_kinetic])
delta_e = np.empty([n_x,n_kinetic])
delta_t = np.empty([n_x,n_kinetic])

# [n_x,n_kinetic,n_time]
for i in range(n_kinetic):
    f   = np.array(sim.moment_zero[:,i,0,:])
    delta_n[:,i] = average(f,t,window,n_x)
    f   = np.array(sim.moment_zero[:,i,1,:])
    delta_e[:,i] = average(f,t,window,n_x)
    delta_t[:,i] = (delta_e[:,i]-1.5*sim.profile['tem_s'][i]*delta_n[:,i])/(1.5*sim.profile['den_s'][i])

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i

color = ['k','m','b','c']

#======================================
fig = plt.figure(figsize=(7*3,6))

f = np.empty([n_x,n_kinetic])

for j in range(3):
    ax = fig.add_subplot(1,3,j+1)
    ax.grid(which="majorminor",ls=":")
    ax.grid(which="major",ls=":")
    ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)
    ax.set_title(r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$',fontsize=GFONTSIZE)

    if j == 0:
        f = delta_n
        ax.set_ylabel('dn',color='k',fontsize=GFONTSIZE)
    if j == 1:
        f = delta_e
        ax.set_ylabel('de',color='k',fontsize=GFONTSIZE)
    if j == 2:
        f = delta_t
        ax.set_ylabel('dt',color='k',fontsize=GFONTSIZE)

    for i in range(n_kinetic):
        stag = sim.tagspec[i]
        ax.plot(sim.profile['r'],f[:,i],label=stag,color=color[i])

    ax.legend()

if ftype == 'screen':
    plt.show()
else:
    outfile = 'moment_zero.'+ftype
    plt.savefig(outfile)
