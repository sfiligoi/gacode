import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from gyro.data_plot import gyrodata_plot

data_in = gyrodata_plot('./')

# Use first 3 args to define plot and font size 
rc('text',usetex=True)
rc('font',size=int(sys.argv[1]))
data_in.lx = int(sys.argv[2])
data_in.ly = int(sys.argv[3])

# Shift list by 3
sys.argv = sys.argv[3:]

plot_type = sys.argv[1]

if plot_type == 'freq':

   w      = float(sys.argv[2])
   ext  = sys.argv[3]

   data_in.plot_freq(w=w)

   outfile = 'freq.'+ext

elif plot_type == 'ball':

   index  = sys.argv[2]
   ext  = sys.argv[3]
   tmax   = float(sys.argv[4])

   key = data_in.plot_balloon(index=index,tmax=tmax)

   outfile = key+'.'+ext

elif plot_type == 'zf':
   
   w      = float(sys.argv[2])
   ext  = sys.argv[3]
   
   data_in.plot_zf(w=w)

   outfile = 'zf.'+ext

elif plot_type == 'flux':

   w        = float(sys.argv[2])
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   ymin     = sys.argv[5]
   ymax     = sys.argv[6]
   ext      = sys.argv[7]
   
   data_in.plot_gbflux(w=w,field=field,i_moment=i_moment,ymin=ymin,ymax=ymax)

   outfile = 'flux.'+ext

elif plot_type == 'flux_i':

   w        = float(sys.argv[2])
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   ymin     = sys.argv[5]
   ymax     = sys.argv[6]
   ext      = sys.argv[7]
   
   data_in.plot_gbflux_i(w=w,field=field,i_moment=i_moment,ymin=ymin,ymax=ymax)

   outfile = 'flux_i.'+ext

elif plot_type == 'flux_n':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   datafile = sys.argv[6]
   ext      = sys.argv[7]
   
   data_in.plot_gbflux_n(field=field,i_moment=i_moment,w=w,datafile=datafile)

   outfile = 'flux_n.'+ext

elif plot_type == 'flux_rt':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   ext      = sys.argv[6]
   
   data_in.plot_flux_rt(field=field,i_moment=i_moment,w=w)

   outfile = 'flux_rt.'+ext

elif plot_type == 'flux_exc':

   simdir   = sys.argv[2]
   w        = float(sys.argv[3])
   ext      = sys.argv[4]
   
   data_in.plot_flux_exc(w=w)

   outfile = 'flux_exc.'+ext

elif plot_type == 'phi_n0':

   simdir = sys.argv[2]
   lx     = float(sys.argv[3])
   ly     = float(sys.argv[4])
   ymax   = sys.argv[5]
   span1  = float(sys.argv[6])
   span2  = float(sys.argv[7])
   ext    = sys.argv[8]

   ax = data_in.plot_phi_n0(lx=lx,ly=ly,ymax=ymax,span1=span1,span2=span2)
   
   outfile = 'phi_n0.'+ext
  
elif plot_type == 'hmom':

   w         = float(sys.argv[2])
   i_kinetic = int(sys.argv[3])
   i_moment  = int(sys.argv[4])
   ext       = sys.argv[5]

   ax = data_in.plot_moment_zero(w=w,i_kinetic=i_kinetic,i_moment=i_moment)
   
   outfile = 'moment_zero.'+ext
  
#---------------------------------------------------------------
# Plot to screen or to image file
if ext == 'screen':
    plt.show()
else:
    plt.savefig(outfile)

if ext != 'screen':
       print 'INFO: (data_plot_single) Created '+outfile
#---------------------------------------------------------------
