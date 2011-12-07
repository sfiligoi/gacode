"""Generate a plot of a profile function"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.profiles_gen.data import profiles_genData

GFONTSIZE=16

prof = profiles_genData(sys.argv[1])
keys = sorted(prof.data.keys())

tag  = sys.argv[2]
rvar = sys.argv[3]

success = 0
for i in range(len(keys)):
    if tag == keys[i].split()[0]:
        success = 1
        fulltag = keys[i]

if success == 0:
    print "ERROR (profiles_gen_plot): Bad profile = "+tag
    sys.exit()

fig = plt.figure(figsize=(12,7))
ax  = fig.add_subplot(111)

if rvar == "r":
    ax.set_xlabel(r"$r \, (m)$",fontsize=GFONTSIZE)
    x = prof.data['rmin (m)']

if rvar == "r/a":
    ax.set_xlabel(r"$r/a$",fontsize=GFONTSIZE)
    x = prof.data['rmin (m)']
    x = x/max(x)

if rvar == "rho":
    ax.set_xlabel(r"$\rho$",fontsize=GFONTSIZE)
    x = prof.data['rho (-)']

ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_ylabel(tag,color='k',fontsize=GFONTSIZE)
ax.plot(x,prof.data[fulltag])
plt.show()
#=====================================
