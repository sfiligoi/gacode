import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from cgyro.data_plot import cgyrodata_plot

rc('text',usetex=True)
rc('font',size=18)

plot_type = sys.argv[1]

if plot_type == 'freq':
   
   w     = float(sys.argv[2])
   ftype = sys.argv[3]

   cgyrodata_plot('./').plot_freq(w=w)

   outfile = 'out.cgyro.freq.'+ftype

elif plot_type == 'ky_freq':
   
   w     = float(sys.argv[2])
   ftype = sys.argv[3]

   cgyrodata_plot('./').plot_ky_freq(w=w)

   outfile = 'out.cgyro.ky_freq.'+ftype

elif plot_type == 'ky_phi':
   
   field = int(sys.argv[2])
   ymin  = sys.argv[3]
   ymax  = sys.argv[4]
   nstr  = sys.argv[5]
   ftype = sys.argv[6]

   cgyrodata_plot('./').plot_ky_phi(field=field,ymin=ymin,ymax=ymax,nstr=nstr)

   outfile = 'out.cgyro.ky_phi.'+ftype

elif plot_type == 'rcorr_phi':
   
   w     = float(sys.argv[2])
   ftype = sys.argv[3]

   cgyrodata_plot('./').plot_rcorr_phi(w=w)

   outfile = 'out.cgyro.rcorr_phi.'+ftype

elif plot_type == 'geo':
   
   ftype = sys.argv[2]

   cgyrodata_plot('./').plot_geo()

   outfile = 'out.cgyro.geo.'+ftype
   
#---------------------------------------------------------------
# Plot to screen or to image file
if ftype == 'screen':
    plt.show()
else:
    print 'INFO: (data_plot_single) Created '+outfile
    plt.savefig(outfile)
#---------------------------------------------------------------
