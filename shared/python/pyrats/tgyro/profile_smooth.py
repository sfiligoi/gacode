import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.tgyro.data import TGYROData
         
sim   = TGYROData(sys.argv[1])
ftype = sys.argv[2]

GFONTSIZE=18
#======================================
fig = plt.figure(figsize=(18,6))
fig.suptitle(sys.argv[1])
fig.subplots_adjust(left=0.05,right=0.95)
fig.subplots_adjust(wspace=0.3)
#=====================================

n   = sim.n_iterations
n_r = sim.n_radial

#----------------------------------------------------------------
# Interpolation onto 32-point grid
#
nf = int(32/n_r)
rf = np.zeros((n_r-1)*nf+1)
zf = np.zeros(((n_r-1)*nf+1,6))
pf = np.zeros(((n_r-1)*nf+1,6))

j = 0
for i in range(n_r-1):
    for m in range(nf):
        z = m/(1.0*nf)
        rf[j] = sim.data['r/a'][n][i]*(1-z)+sim.data['r/a'][n][i+1]*z
        zf[j,0] = sim.data['a/Lne'][0][i]*(1-z)+sim.data['a/Lne'][0][i+1]*z
        zf[j,1] = sim.data['a/LTi'][0][i]*(1-z)+sim.data['a/LTi'][0][i+1]*z
        zf[j,2] = sim.data['a/LTe'][0][i]*(1-z)+sim.data['a/LTe'][0][i+1]*z
        zf[j,3] = sim.data['a/Lne'][n][i]*(1-z)+sim.data['a/Lne'][n][i+1]*z
        zf[j,4] = sim.data['a/LTi'][n][i]*(1-z)+sim.data['a/LTi'][n][i+1]*z
        zf[j,5] = sim.data['a/LTe'][n][i]*(1-z)+sim.data['a/LTe'][n][i+1]*z
        j = j+1

rf[j] = sim.data['r/a'][n][n_r-1]

zf[j,0] = sim.data['a/Lne'][0][n_r-1]
zf[j,1] = sim.data['a/LTi'][0][n_r-1]
zf[j,2] = sim.data['a/LTe'][0][n_r-1]
zf[j,3] = sim.data['a/Lne'][n][n_r-1]
zf[j,4] = sim.data['a/LTi'][n][n_r-1]
zf[j,5] = sim.data['a/LTe'][n][n_r-1]

pf[j,0] = sim.data['ne'][0][n_r-1]
pf[j,1] = sim.data['ti'][0][n_r-1]
pf[j,2] = sim.data['te'][0][n_r-1]
pf[j,3] = sim.data['ne'][n][n_r-1]
pf[j,4] = sim.data['ti'][n][n_r-1]
pf[j,5] = sim.data['te'][n][n_r-1]

# Exponential integration to obtain smooth profiles
for i in np.arange(j,0,-1):
    pf[i-1,:] = pf[i,:]*np.exp(0.5*(rf[i]-rf[i-1])*(zf[i,:]+zf[i-1,:]))

#----------------------------------------------------------------

ax = fig.add_subplot(131)
ax.plot(rf,pf[:,0],label='init')
ax.plot(rf,pf[:,3],label='final')
ax.set_ylabel('$n_e$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(132)
ax.plot(rf,pf[:,1],label='init')
ax.plot(rf,pf[:,4],label='final')
ax.set_ylabel('$T_i$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(133)
ax.plot(rf,pf[:,2],label='init')
ax.plot(rf,pf[:,5],label='final')
ax.set_ylabel('$T_e$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()


if ftype == 'screen':
    plt.show()
else:
    outfile = 'profile.'+ftype
    plt.savefig(outfile)
