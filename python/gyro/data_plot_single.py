import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data_plot import gyrodata_plot

plot_type = sys.argv[1]

if plot_type == 'freq':

   simdir = sys.argv[2]
   w      = float(sys.argv[3])
   ftype  = sys.argv[4]

   gyrodata_plot(simdir).plot_freq(w=w)

   outfile = 'freq.'+ftype

elif plot_type == 'balloon':

   simdir = sys.argv[2]
   index  = sys.argv[3]
   ftype  = sys.argv[4]
   tmax   = float(sys.argv[5])

   key = gyrodata_plot(simdir).plot_balloon(index=index,tmax=tmax)

   outfile = key+'.'+ftype

elif plot_type == 'zf':

   simdir = sys.argv[2]
   w      = float(sys.argv[3])
   ftype  = sys.argv[4]
   
   gyrodata_plot(simdir).plot_zf(w=w)

   outfile = 'zf.'+ftype

elif plot_type == 'gbflux':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   lx       = float(sys.argv[6])
   ly       = float(sys.argv[7])
   title    = sys.argv[8]
   ymin     = sys.argv[9]
   ymax     = sys.argv[10]
   span1   = float(sys.argv[11])
   span2   = float(sys.argv[12])
   ftype    = sys.argv[13]
   
   gyrodata_plot(simdir).plot_gbflux(field=field,i_moment=i_moment,w=w,
                                     lx=lx,ly=ly,title=title,ymin=ymin,ymax=ymax,
                                     span1=span1)

   outfile = 'gbflux.'+ftype

elif plot_type == 'phi_n0':

   simdir = sys.argv[2]
   lx     = float(sys.argv[3])
   ly     = float(sys.argv[4])
   ymax   = sys.argv[5]
   span1  = float(sys.argv[6])
   span2  = float(sys.argv[7])
   ftype  = sys.argv[8]

   ax = gyrodata_plot(simdir).plot_phi_n0(lx=lx,ly=ly,ymax=ymax,span1=span1,span2=span2)
   
   outfile = 'phi_n0.'+ftype
  

#---------------------------------------------------------------
# Plot to screen or to image file
if ftype == 'screen':
    plt.show()
else:
    plt.savefig(outfile)
#---------------------------------------------------------------
