import sys
import numpy as np
import matplotlib.pyplot as plt
from tgyro.data import TGYROData
         
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
ax.plot(sim.data['r/a'][0],sim.data['ne'][0],label='init')
ax.plot(sim.data['r/a'][0],sim.data['ne'][n],label='final')
ax.set_ylabel('$n_e$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(132)
ax.plot(sim.data['r/a'][0],sim.data['ti'][0],label='init')
ax.plot(sim.data['r/a'][0],sim.data['ti'][n],label='final')
ax.set_ylabel('$T_i$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

ax = fig.add_subplot(133)
ax.plot(sim.data['r/a'][0],sim.data['te'][0],label='init')
ax.plot(sim.data['r/a'][0],sim.data['te'][n],label='final')
ax.set_ylabel('$T_e$',color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.legend()

if ftype == 'screen':
    plt.show()
else:
    outfile = 'profile.'+ftype
    plt.savefig(outfile)
