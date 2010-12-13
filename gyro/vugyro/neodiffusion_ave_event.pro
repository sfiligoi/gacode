pro neodiffusion_ave_event, neodiffusion_ave

  common GLOBAL
  common TAG_DATA
  common PLOT_VARIABLES
  common PRIVATE_NEODIFFUSION_AVE
  
  widget_control, neodiffusion_ave.id, $
                  get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: begin
        i_spec = i_spec+1
        if (i_spec gt n_kinetic-1) then i_spec = 0
        goto, plot_it
     end   

     2: begin
        i_neomom = i_neomom+1
        if (i_neomom gt 2) then i_neomom = 0
        goto, plot_it
     end   

     3: begin
        i_runave = 1-i_runave
        goto, plot_it
     end

     4: begin
        zoom = zoom*2
        goto, plot_it
     end

     5: begin
        zoom = zoom/2
        goto, plot_it
     end

     6: begin
        i_units = 1-i_units
        goto, plot_it
     end

     7: begin
        plot_mode = 2
        goto, plot_it
     end

     8: widget_control, neodiffusion_ave.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  y = fltarr(n_time)
  y_tave = fltarr(n_time)

  tag_gen,n_ion,i_spec,i_neomom,0,0

  if i_neomom eq 2 then begin
     title = 'parallel velocity/cs'
     pname = 'neodiff-'+ftag_s+'_vpar'
   endif else begin
     title = tag_s+' '+tag_m+' Diffusion_neo'
     pname = 'neodiff-'+ftag_s+ftag_m
  endelse

  ytitle = title

  y[*] = diff_neo[i_spec,i_neomom,*]

  if (i_units eq 0) then begin
     plot_units = 1.0
     units = ' [units of (c!ds!n/a)!4q!3!ds!n!u2!n]'
  endif else begin
     plot_units = xunits[8]
     units = ' [units of (m!u2!n/sec)]'
  endelse

  plot_def_new,pname

  ;; Set units here
  y[*] = plot_units*y[*]

  y_axis_magic,y,ymin,ymax,d_y

  diff_statistics,y,t_min,t_max,i1,i2,A_f,D_RMS

  unity = fltarr(i2-i1)
  unity[*] = 1.0

  A_r = strmid(strtrim(string(A_f),2),0,5)
  D_r = strmid(strtrim(string(D_RMS),2),0,5)

  text = ' !3['+A_r+'!9+!3'+D_r+']'

  if (i_runave eq 1) then begin

     y_tave[0] = 0.0
     for i_t = 1,n_time-1 do begin
        t_total = t[i_t]-t[0]
        sum_f  = 0.0
        for i=0,i_t-1 do begin
           dt = t[i+1]-t[i]
           f  = 0.5*(y[i+1]+y[i])
           sum_f  = sum_f + f*dt
        endfor
        y_tave[i_t] = sum_f/t_total
     endfor

  endif

  plot,[0],[0],$
       /nodata,$
       title=title+text,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(t),max(t)],$
       xtitle='!3(c!ds!n/a) t',$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax]/zoom,$
       ytitle=ytitle+units, $
       color=line

  ;; Diffusion trace

  if (i_runave eq 0) then begin
     oplot,t,y,color=color_vec[0]
  endif else begin
     oplot,t,y_tave,color=color_vec[0]
  endelse
  
  ;; Average line

  oplot,t[i1:i2],unity*A_f,color=color_vec[1]   

  ;; Curve dots

  oplot,t[i1]*[1,1],y[i1]*[1,1],psym=8,color=line   
  oplot,t[i2]*[1,1],y[i2]*[1,1],psym=8,color=line   

  ;; RMS deviation bars

  oplot,t[i1:i2],unity*(A_f-D_RMS),color=color_vec[2]   
  oplot,t[i1:i2],unity*(A_f+D_RMS),color=color_vec[2]   

  plot_finish

  return

end
