import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data_plot import gyrodata_plot

fig = plt.figure(figsize=(16,8))

GFONTSIZE=16

sim = gyrodata_plot(sys.argv[1])
sim.plot_freq(window=float(sys.argv[2]),fig=fig)
ftype = sys.argv[3]

if ftype == 'screen':
    plt.show()
else:
    outfile = 'freq.'+ftype
    plt.savefig(outfile)
