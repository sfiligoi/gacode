import sys
from gacodeplotdefs import *
from gyro.data_plot import gyrodata_plot

sim    = gyrodata_plot(sys.argv[1])
key = sim.plot_balloon(index=sys.argv[2], tmax=float(sys.argv[5]))
ftype  = sys.argv[3]

if ftype == 'screen':
    plt.show()
else:
    outfile = key+'.'+ftype
    plt.savefig(outfile)
