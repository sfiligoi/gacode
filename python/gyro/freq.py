import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data_plot import gyrodata_plot

# Command-line args
simdir = sys.argv[1]
w      = float(sys.argv[2])
ftype  = sys.argv[3]

fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.09,right=0.97,top=0.92,bottom=0.12)

gyrodata_plot(simdir).plot_freq(window=w,fig=fig)

if ftype == 'screen':
    plt.show()
else:
    outfile = 'freq.'+ftype
    plt.savefig(outfile)
