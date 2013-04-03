import sys
import numpy as np
from gacodeplotdefs import *
from pyrats.tgyro.data import TGYROData

sim   = TGYROData(sys.argv[1])
ftype = sys.argv[2]
dots  = int(sys.argv[3])

#======================================
fig = plt.figure(figsize=(12,10))
fig.suptitle('\\verb|'+sys.argv[1]+'|')
fig.subplots_adjust(top=0.93,bottom=0.07)
fig.subplots_adjust(left=0.07,right=0.96)
fig.subplots_adjust(wspace=0.3)
#=====================================

n   = sim.n_iterations
n_r = sim.n_radial

#----------------------------------------------------------------
# Interpolation onto 32-point grid
#
nf = int(64/n_r)
rf = np.zeros((n_r-1)*nf+1)
zf = np.zeros(((n_r-1)*nf+1,8))
pf = np.zeros(((n_r-1)*nf+1,8))

a = sim.data['rmin'][0][1]/sim.data['r/a'][0][1]

#----------------------------------------------------------------------------
# Want to plot Omega, the rotation frequency.  This is cumbersome.
#
# Note conversion of c_s, below, from m/s to cm/s
#
# Derivative on sparse TGYRO grid: a*d(w0)/dr
aw0p = -sim.data['a*gamma_p/cs']*(100*sim.data['c_s'])/a/sim.data['rmaj/a']
# Rotation on sparse TGYRO grid: w0
w0 = sim.data['M=wR/cs']*(100*sim.data['c_s'])/(a*sim.data['rmaj/a'])

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
        zf[j,6] = aw0p[0][i]*(1-z)+aw0p[0][i+1]*z
        zf[j,7] = aw0p[n][i]*(1-z)+aw0p[n][i+1]*z
        j = j+1

rf[j] = sim.data['r/a'][n][n_r-1]

zf[j,0] = sim.data['a/Lne'][0][n_r-1]
zf[j,1] = sim.data['a/LTi'][0][n_r-1]
zf[j,2] = sim.data['a/LTe'][0][n_r-1]
zf[j,3] = sim.data['a/Lne'][n][n_r-1]
zf[j,4] = sim.data['a/LTi'][n][n_r-1]
zf[j,5] = sim.data['a/LTe'][n][n_r-1]
zf[j,6] = aw0p[0][n_r-1]
zf[j,7] = aw0p[n][n_r-1]

pf[j,0] = sim.data['ne'][0][n_r-1]
pf[j,1] = sim.data['ti'][0][n_r-1]
pf[j,2] = sim.data['te'][0][n_r-1]
pf[j,3] = sim.data['ne'][n][n_r-1]
pf[j,4] = sim.data['ti'][n][n_r-1]
pf[j,5] = sim.data['te'][n][n_r-1]
pf[j,6] = w0[0][n_r-1]
pf[j,7] = w0[n][n_r-1]


# Exponential integration to obtain smooth profiles
for i in np.arange(j,0,-1):
    pf[i-1,0:6] = pf[i,0:6]*np.exp(0.5*(rf[i]-rf[i-1])*(zf[i,0:6]+zf[i-1,0:6]))
# Linear integration for rotation
for i in np.arange(j,0,-1):
    pf[i-1,6:8] = pf[i,6:8]-0.5*(rf[i]-rf[i-1])*(zf[i,6:8]+zf[i-1,6:8])

#----------------------------------------------------------------

ax = fig.add_subplot(221)
ax.plot(rf,pf[:,0],label='init')
ax.plot(rf,pf[:,3],label='final')
if dots == 1:
    ax.plot(sim.data['r/a'][n],sim.data['ne'][n],'o',color='k')
ax.set_ylabel(r'$n_e$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(222)
ax.plot(rf,pf[:,1],label='init')
ax.plot(rf,pf[:,4],label='final')
if dots == 1:
    ax.plot(sim.data['r/a'][n],sim.data['ti'][n],'o',color='k')
ax.set_ylabel(r'$T_i$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(223)
ax.plot(rf,pf[:,2],label='init')
ax.plot(rf,pf[:,5],label='final')
if dots == 1:
    ax.plot(sim.data['r/a'][n],sim.data['te'][n],'o',color='k')
ax.set_ylabel(r'$T_e$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(224)
ax.plot(rf,pf[:,6]/1e3,label='init')
ax.plot(rf,pf[:,7]/1e3,label='final')
if dots == 1:
    ax.plot(sim.data['r/a'][n],w0[n]/1e3,'o',color='k')
ax.set_ylabel(r'$\omega_0\;\mathrm{(krad/s)}$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)

ax.legend()


if ftype == 'screen':
    plt.show()
else:
    outfile = 'profile.'+ftype
    plt.savefig(outfile)
