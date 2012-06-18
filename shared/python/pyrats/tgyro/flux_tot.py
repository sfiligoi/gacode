import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.tgyro.data import TGYROData
         
sim   = TGYROData(sys.argv[1])
units = int(sys.argv[2])
ftype = sys.argv[3]

GFONTSIZE=18
#======================================
fig = plt.figure(figsize=(18,6))
fig.suptitle(sys.argv[1])
fig.subplots_adjust(left=0.05,right=0.95)
fig.subplots_adjust(wspace=0.3)
#=====================================

n = sim.n_iterations

x   = sim.data['r/a'][0]
ggb = sim.data['Gamma_GB'][n]
qgb = sim.data['Q_GB'][n]

ax = fig.add_subplot(131)
if units == 0:
    ax.plot(x,sim.data['pflux_e_tot'][n],label='total')
    ax.plot(x,sim.data['pflux_e_target'][n],label='target')
    ax.set_ylabel('$\Gamma/\Gamma_{GB}$',color='k',fontsize=GFONTSIZE)
else:
    ax.plot(x,sim.data['pflux_e_tot'][n]*ggb,label='total')
    ax.plot(x,sim.data['pflux_e_target'][n]*ggb,label='target')
    ax.set_ylabel('$\Gamma [10^{19}/m^2/s]$',color='k',fontsize=GFONTSIZE)
    
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(132)
if units == 0:
    ax.plot(x,sim.data['eflux_e_tot'][n],label='total')
    ax.plot(x,sim.data['eflux_e_target'][n],label='target')
    ax.set_ylabel('$Q_e/Q_{GB}$',color='k',fontsize=GFONTSIZE)
else:
    ax.plot(x,sim.data['eflux_e_tot'][n]*qgb,label='total')
    ax.plot(x,sim.data['eflux_e_target'][n]*qgb,label='target')
    ax.set_ylabel('$Q_e [MW/m^2]$',color='k',fontsize=GFONTSIZE)

ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(133)
if units == 0:
    ax.plot(x,sim.data['eflux_i_tot'][n],label='total')
    ax.plot(x,sim.data['eflux_i_target'][n],label='target')
    ax.set_ylabel('$Q_i/Q_{GB}$',color='k',fontsize=GFONTSIZE)
else:
    ax.plot(x,sim.data['eflux_i_tot'][n]*qgb,label='total')
    ax.plot(x,sim.data['eflux_i_target'][n]*qgb,label='target')
    ax.set_ylabel('$Q_i [MW/m^2]$',color='k',fontsize=GFONTSIZE)

ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

if ftype == 'screen':
    plt.show()
else:
    outfile = 'flux_tot.'+ftype
    plt.savefig(outfile)
