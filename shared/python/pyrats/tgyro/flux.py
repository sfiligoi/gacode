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

n = sim.n_iterations

ax = fig.add_subplot(131)
ax.plot(sim.data['r/a'][4],sim.data['pflux_e_tur'][n],label='tur')
ax.plot(sim.data['r/a'][4],sim.data['pflux_e_neo'][n],label='neo')
ax.set_ylabel('$\Gamma/\Gamma_{GB}$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(132)
ax.plot(sim.data['r/a'][4],sim.data['eflux_e_tur'][n],label='tur')
ax.plot(sim.data['r/a'][4],sim.data['eflux_e_neo'][n],label='neo')
ax.set_ylabel('$Q_e/Q_{GB}$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(133)
ax.plot(sim.data['r/a'][4],sim.data['eflux_i_tur'][n],label='tur')
ax.plot(sim.data['r/a'][4],sim.data['eflux_i_neo'][n],label='neo')
ax.set_ylabel('$Q_i/Q_{GB}$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

if ftype == 'screen':
    plt.show()
else:
    outfile = 'flux.'+ftype
    plt.savefig(outfile)
