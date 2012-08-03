import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.tgyro.data import TGYROData
         
sim   = TGYROData(sys.argv[1])
f     = sys.argv[2]
ftype = sys.argv[3]

GFONTSIZE=18
#======================================
fig = plt.figure(figsize=(12,10))
fig.suptitle(sys.argv[1])
#=====================================

n = sim.n_iterations

ax = fig.add_subplot(111)
ax.set_ylabel(f,color='k',fontsize=GFONTSIZE)
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)

for i in range(n+1):
    ax.plot(sim.data['r/a'][i],sim.data[f][i],label=str(i))

ax.legend()

if ftype == 'screen':
    plt.show()
else:
    outfile = 'profile.'+ftype
    plt.savefig(outfile)
