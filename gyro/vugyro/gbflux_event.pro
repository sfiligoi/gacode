pro gbflux_event, gbflux_wid

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_GBFLUX  

  widget_control, gbflux_wid.id, get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     0: goto, plot_it

     1: begin
        i_f = i_f+1
        if (i_f ge n_field) then i_f = -1
        goto, plot_it
     end

     2: begin
        i_spec = i_spec+1
        if (i_spec ge n_kinetic) then i_spec = 0
        goto, plot_it
     end   

     3: begin
        i_moment = i_moment+1
        if (i_moment ge p_moment) then i_moment = 0
        goto, plot_it
     end   

     4: begin
        zoom = 2*zoom
        goto, plot_it
     end

     5: begin
        zoom = zoom/2
        goto, plot_it
     end

     6: begin
        if (exists_exp_derived eq 1) then begin
           i_units = i_units+1
           if (i_units eq 3) then i_units = 0
        endif
        goto, plot_it
     end

     7: begin
        if (0 eq 1) then begin
           i_tp = i_tp+1
           if i_tp ge 3 then i_tp = 0
        endif
        goto, plot_it
     end

     8: begin
        plot_mode = 2
        goto, plot_it
     end

     9: begin
        plot_export = 1
        goto, plot_it
     end

     10: widget_control, gbflux_wid.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  y = fltarr(n_time)
  
  gbflux_tag,title,pname,ytitle,units1,units2,units3

  case (i_units) of

  0: begin
     ;; flux in gyrobohm units
     xnorm = 1.0
     units = units1
     end

  1: begin
     ;; flux in W/m^2, etc.
     xnorm = xunits[9+i_moment]
     units = units2
     end 

  2: begin 
     ;; flux in W, etc.
     xnorm = xunits[9+i_moment]
     units = units3
     end

  endcase

  case (i_tp) of

     0: begin

        pname = 'gbflux-'+pname
        title = title
        if i_f ge 0 then begin
           y[*] = gbflux[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = 0.0
           for ix=0,n_field-1 do begin
              y[*] = y[*]+gbflux[i_spec,ix,i_moment,*]
           endfor
        endelse

     end

     1: begin

        pname = 'gbflux_t-'+pname
        title = title+' (trapped)'
        if i_f ge 0 then begin
           y[*] = gbflux_trapped[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = 0.0
           for ix=0,n_field-1 do begin 
              y[*] = y[*]+gbflux_trapped[i_spec,ix,i_moment,*]
           endfor
        endelse

     end

     2: begin

        pname = 'gbflux_p-'+pname
        title = title+' (passing)'
        if i_f ge 0 then begin
           y[*] = gbflux[i_spec,i_f,i_moment,*]-$
                  gbflux_trapped[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = 0.0
           for ix=0,n_field-1 do begin
              y[*] = y[*] + gbflux[i_spec,ix,i_moment,*]-$
                     gbflux_trapped[i_spec,ix,i_moment,*] 
           endfor
        endelse

     end

  endcase

  plot_def_new,pname

  ;; Set units here

  y[*] = y[*]*xnorm
  if (i_units eq 2) then begin
     ;; Multiply by V'
     y[*] = y[*]*INTERPOL(exp_derived[22,*],r_from_rho,r[n_r/2])
  endif

  y_axis_magic,y,ymin,ymax,d_y

  unity = fltarr(it2-it1)
  unity[*] = 1.0

  if (it2-it1) gt 20 then begin

     ;; Enough points for Mikkelsen's method

     res = meanstddev(pname,0,y,t,t_min,t_max)
     ave_y   = res.avg
     err_y   = res.sig2
     trend_y = res.trend
     hump_y  = res.hump

  endif else begin

     ;; Only compute simple average

     diff_stat_fast,y,it1,it2,ave_y
     err_y   = -1
     trend_y = -1
     hump_y  = -1

  endelse

  ave_str   = strmid(strtrim(string(ave_y),2),0,5)
  err_str   = strmid(strtrim(string(err_y),2),0,5)
  trend_str = strmid(strtrim(string(trend_y),2),0,5)
  hump_str  = strmid(strtrim(string(hump_y),2),0,5)

  text = ' !3['+ave_str+'!9+!3'+err_str+']'

  plot,[0],[0],$
       /nodata,$
       title=title+text,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(t),max(t)],$
       xtitle=csa_string,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax]/zoom,$
       ytitle=ytitle+units, $
       color=line
  
  ;; flux trace

  oplot,t,y,color=color_vec[0]

  ;; Average line

  oplot,t[it1:it2],unity*ave_y,color=color_vec[1]   

  ;; Curve dots

  oplot,t[it1]*[1,1],y[it1]*[1,1],psym=8,color=line   
  oplot,t[it2]*[1,1],y[it2]*[1,1],psym=8,color=line   

  ;; RMS deviation bars

  oplot,t[it1:it2],unity*(ave_y-err_y),color=color_vec[2]   
  oplot,t[it1:it2],unity*(ave_y+err_y),color=color_vec[2]   

  dy = ymax-ymin

  xyouts,0.7*max(t),(ymin+0.95*dy)/zoom,$
         '!3trend='+trend_str

  xyouts,0.7*max(t),(ymin+0.9*dy)/zoom,$
         '!3hump='+hump_str

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     printf,1,ave_y,err_y
     printf,1,' Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------
  
  plot_finish

  return

end
