"""Generate a plot of a profile function"""

import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from pyrats.profiles_gen.data import profiles_genData

GFONTSIZE=16

rvar    = sys.argv[1]
infiles = sys.argv[2]
plots   = sys.argv[3]

plotvec = string.splitfields(plots,',')

n = len(plotvec)

tag  = sys.argv[2]

for j in range(n):

    print infiles[j]
    prof = profiles_genData(infiles[j])
    keys = sorted(prof.data.keys())

    success = 0
    for i in range(len(keys)):
        if tag == keys[i].split()[0]:
            success = 1
            fulltag = keys[i]

        if success == 0:
            print "ERROR (profiles_gen_plot): Bad profile = "+tag
            sys.exit()

        if j==0:
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
