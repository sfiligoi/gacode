import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data import GYROData

sim       = GYROData(sys.argv[1])
plotfile  = sys.argv[2]
datafile  = sys.argv[3]
lx        = float(sys.argv[4])
ly        = float(sys.argv[5])
title     = sys.argv[6]
ymax      = float(sys.argv[7])

t = sim.t['(c_s/a)t']

# Read data in gbflux_i and make gbflux
sim.read_field_rms()

#======================================
fig = plt.figure(figsize=(lx,ly))
fig.subplots_adjust(left=0.1,right=0.96,top=0.93,bottom=0.13)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$\langle e \phi/T_e \rangle/\rho_\star $',color='k',fontsize=GFONTSIZE)
#=====================================

# rho_* for normalization
x=np.average(sim.profile['rho_s'])

if datafile == 'none':
    # Plot data to screen or image file.
    ax.plot(t,sim.field_rms[0,:]/x,color='k',label=r'$n=0$')
    ax.plot(t,sim.field_rms[1,:]/x,color='purple',label=r'$n>0$')
else:
    # Write data to datafile
    print 'INFO: (gyro_plot) Output to datafile not supported.  Use raw out.gyro.gbflux.'

ax.set_xlim([0,t[-1]])

if ymax > 0:
    ax.set_ylim([0,ymax])

ax.legend()

if plotfile == 'screen':
    plt.show()
else:
    plt.savefig(plotfile)
    print "INFO: (gyro_plot) Wrote plot to " + plotfile + "."
