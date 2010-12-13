pro plot_finish

  common GLOBAL
  common PLOT_VARIABLES

  if (ps_label eq 0) then begin

     xyouts,1,1,/device, $
      '!3'+simdir+' ['+version_tag+'] ['+version_date+']', $
      color=line,$
      size=annosize

  endif

  if (plot_mode eq 2) then begin

     print,'Wrote '+filename
     plot_mode = 1
     device_reset     

     ; reset thicknesses

     !p.charthick=1.0
     !p.thick=1.0
     !x.thick=1.0
     !y.thick=1.0

  endif

end

