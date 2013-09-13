
         integer  
     .         nbeam_pts,nplt_max
c
      parameter (nplt_max = 3000 )
      real *8
     .      timplot(nplt_max),waveform(nplt_max),
     .      s0_start,scale_factor,eps
      common /plot_beam/
     .      waveform,timplot,s0_start,
     .      scale_factor,nbeam_pts
