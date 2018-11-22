import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from gyro.data_plot import gyrodata_plot

data_in = gyrodata_plot('./')

rc('text',usetex=True)
rc('font',size=18)

plot_type = sys.argv[1]

if plot_type == 'freq':

   w      = float(sys.argv[2])
   ftype  = sys.argv[3]

   data_in.plot_freq(w=w)

   outfile = 'freq.'+ftype

elif plot_type == 'balloon':

   index  = sys.argv[2]
   ftype  = sys.argv[3]
   tmax   = float(sys.argv[4])

   key = data_in.plot_balloon(index=index,tmax=tmax)

   outfile = key+'.'+ftype

elif plot_type == 'zf':

   simdir = sys.argv[2]
   w      = float(sys.argv[3])
   ftype  = sys.argv[4]
   
   data_in.plot_zf(w=w)

   outfile = 'zf.'+ftype

elif plot_type == 'gbflux':

   field    = sys.argv[2]
   i_moment = int(sys.argv[3])
   w        = float(sys.argv[4])
   lx       = float(sys.argv[5])
   ly       = float(sys.argv[6])
   title    = sys.argv[7]
   ymin     = sys.argv[8]
   ymax     = sys.argv[9]
   span1   = float(sys.argv[10])
   span2   = float(sys.argv[11])
   ftype    = sys.argv[12]
   
   data_in.plot_gbflux(field=field,i_moment=i_moment,w=w,lx=lx,ly=ly,
                       title=title,ymin=ymin,ymax=ymax,span1=span1)

   outfile = 'gbflux.'+ftype

elif plot_type == 'gbflux_i':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   ymin     = sys.argv[6]
   ymax     = sys.argv[7]
   ftype    = sys.argv[8]
   
   data_in.plot_gbflux_i(field=field,i_moment=i_moment,w=w,ymin=ymin,ymax=ymax)

   outfile = 'gbflux_i.'+ftype

elif plot_type == 'gbflux_n':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   datafile = sys.argv[6]
   ftype    = sys.argv[7]
   
   data_in.plot_gbflux_n(field=field,i_moment=i_moment,w=w,datafile=datafile)

   outfile = 'gbflux_n.'+ftype

elif plot_type == 'gbflux_rt':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   ftype    = sys.argv[6]
   
   data_in.plot_gbflux_rt(field=field,i_moment=i_moment,w=w)

   outfile = 'gbflux_rt.'+ftype

elif plot_type == 'gbflux_exc':

   simdir   = sys.argv[2]
   w        = float(sys.argv[3])
   ftype    = sys.argv[4]
   
   data_in.plot_gbflux_exc(w=w)

   outfile = 'gbflux_exc.'+ftype

elif plot_type == 'phi_n0':

   simdir = sys.argv[2]
   lx     = float(sys.argv[3])
   ly     = float(sys.argv[4])
   ymax   = sys.argv[5]
   span1  = float(sys.argv[6])
   span2  = float(sys.argv[7])
   ftype  = sys.argv[8]

   ax = data_in.plot_phi_n0(lx=lx,ly=ly,ymax=ymax,span1=span1,span2=span2)
   
   outfile = 'phi_n0.'+ftype
  
elif plot_type == 'moment_zero':

   i_kinetic = int(sys.argv[2])
   w         = float(sys.argv[3])
   ftype     = sys.argv[4]

   ax = data_in.plot_moment_zero(i_kinetic=i_kinetic,w=w)
   
   outfile = 'moment_zero.'+ftype

elif plot_type == 'source':
   
   i_kinetic = int(sys.argv[2])
   w        = float(sys.argv[3])
   ftype    = sys.argv[4]

   ax = data_in.plot_source(i_kinetic=i_kinetic,w=w)
   
   outfile = 'out.gyro.source.'+ftype
  
#---------------------------------------------------------------
# Plot to screen or to image file
if ftype == 'screen':
    plt.show()
else:
    plt.savefig(outfile)

if ftype != 'screen':
       print 'INFO: (data_plot_single) Created '+outfile
#---------------------------------------------------------------
