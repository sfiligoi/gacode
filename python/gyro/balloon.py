import sys
from gacodeplotdefs import *
from gyro.data_plot import gyrodata_plot

# Command-line args
simdir = sys.argv[1]
index  = sys.argv[2]
ftype  = sys.argv[3]
tmax   = float(sys.argv[4])

key = gyrodata_plot(simdir).plot_balloon(index=index,tmax=tmax)

if ftype == 'screen':
    plt.show()
else:
    outfile = key+'.'+ftype
    plt.savefig(outfile)
