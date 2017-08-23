import sys
from gacodeplotdefs import *
from gyro.data_plot import gyrodata_plot

fig = plt.figure(figsize=(10,6))
sim   = gyrodata_plot(sys.argv[1])
w     = float(sys.argv[2])
sim.plot_zf(w=w,fig=fig)

ftype = sys.argv[3]
if ftype == 'screen':
    plt.show()
else:
    outfile = 'zf.'+ftype
    plt.savefig(outfile)

