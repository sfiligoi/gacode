pro plot_def_new,label

  common GLOBAL
  common PLOT_VARIABLES

  ;;----------------------------------
  ;; Some common plotting definitions:
  ;;  
  if (plot_mode eq 1) then begin

     ;;screen;

     line = 1
     
     color_vec = [5,6,7,10,2,3]
     line_vec  = [0,0,0,1,0,0]

     !p.charsize = c_size
     dotsize=0.5
     
  endif else begin

     ;;postscript;

     line = 0

     color_vec = [4,6,8,12,2,3]
     if ps_color eq 0 then color_vec = [0,0,0,0,0,0]

     line_vec  = [0,0,0,1,0,0]
     if ps_color eq 0 then line_vec  = [0,1,0,0,0,0]

     ;;Postscript parameters:;
     
     set_plot,'ps'

     if ps_val eq 0 then begin

        filename = label+'.ps'

        device,/color, $
               bits=8, $
               filename=filename, $
               xsize=xsize, $
               ysize=ysize

     endif else begin

        filename = label+'.eps'

        device,/color, $
               /encapsulated, $
               bits=8, $
               filename=filename, $
               xsize=xsize, $
               ysize=ysize

     endelse

     !p.charsize = 1.25
     dotsize=0.25

     if (ps_thick eq 1) then begin

        ;;thick;

        !p.charthick=3.5
        !p.thick=3.0
        !x.thick=3.0
        !y.thick=3.0

     endif

  endelse
  ;;----------------------------------;

  return

end
