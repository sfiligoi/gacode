pro sdiffusion_ave_event, sdiffusion_ave

  common GLOBAL
  common TAG_DATA
  common PLOT_VARIABLES
  common PRIVATE_SDIFFUSION_AVE  
  
  widget_control, sdiffusion_ave.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: begin
        i_spec = i_spec+1
        if (i_spec ge n_kinetic) then i_spec = 0
        goto, plot_it
     end   

     2: begin
        i_moment = 1+i_moment
        if (i_moment ge 2) then i_moment = 0
        goto, plot_it
     end   

     3: begin
        zoom = zoom*2
        goto, plot_it
     end

     4: begin
        zoom = zoom/2
        goto, plot_it
     end

     5: begin
        i_units = 1-i_units
        goto, plot_it
     end

     6: begin
        if (exists_sdiff_t eq 1) then begin
           i_tp = i_tp+1
           if i_tp ge 3 then i_tp = 0
        endif
        goto, plot_it
     end

     7: begin
        plot_mode = 2
        goto, plot_it
     end

     8: widget_control, sdiffusion_ave.top, /destroy

     20: begin
        i_ptype = 0
        goto, plot_it
     end

     21: begin
        i_ptype = 1
        goto, plot_it
     end

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  y = fltarr(n_time)

  y1 = fltarr(n_time)
  y2 = fltarr(n_time)
  y3 = fltarr(n_time)
  y4 = fltarr(n_time)
  
  tag_gen,n_ion,i_spec,i_moment,0,0

  title = tag_s
  pname = ftag_s
  if (i_moment eq 0) then begin
     title = title+' Mom.'
     pname = pname+'_mom'
  endif else begin
     title = title+' Heating'
     pname = pname+'_heat'
  endelse
  ytitle = pname

  case (i_tp) of

     0: begin

        pname = 'sdiff-'+pname

        y[*] = sdiff[i_spec,i_moment,*]

        y1[*] = spdiff[i_spec,i_moment,0,*]
        y2[*] = spdiff[i_spec,i_moment,1,*]
        y3[*] = spdiff[i_spec,i_moment,2,*]

        ;;   y4[*] =y1[*]+y2[*]+y3[*]

      ;;        y[*] = y4[*]
      ;; shouldn't need to do this?

     end

     1: begin

        pname = 'sdiff_t-'+pname
        title = title+' (trap)'
        y[*] = sdiff_trapped[i_spec,i_moment,*]

     end

     2: begin

        pname = 'sdiff_p-'+pname
        title = title+' (pass)'
        y[*] = sdiff[i_spec,i_moment,*]-$
          sdiff_trapped[i_spec,i_moment,*]

     end

  endcase

  if (i_units eq 0) then begin
     plot_units = 1.0
     units = ' [units of '+chi_string+']'
  endif else begin
     plot_units = xunits[8]
     units = ' [units of (m!u2!n/sec)]'
  endelse

  if i_ptype eq 0 then begin

     ;; Diffusion TIME trace

     plot_def_new,pname

     ;; Set units here

     y[*] = plot_units*y[*]

     y1[*] = plot_units*y1[*]
     y2[*] = plot_units*y2[*]
     y3[*] = plot_units*y3[*]

     y4[*] = plot_units*y4[*]
     y4[*] = 1.01*y4[*]

      ymax = 1000.
      ymin = -1000.

     for it = 0,n_time-1 do begin
      if y(it) gt ymax then begin
        y(it) = 1.
      endif
      if y(it) lt ymin then begin
        y(it) = -1.
      endif

      if y(it) eq 0.0 then begin
       y(it) = 0.00001
      endif
     endfor

     y_axis_magic,y,ymin,ymax,d_y

     print , 'ymin=',ymin,'ymax=',ymax
    
    if ymax gt 0. then begin
     ymin = -ymax
    endif
    if ymin lt 0. then begin
     ymax = -ymin
    endif

     unity = fltarr(it2-it1)
     unity[*] = 1.0
  ;;       res = meanstddev(pname,0,y,t,t_min,t_max)
  ;;    ave_y = res.avg
  ;;   err_y = res.sig2

       diff_statistics,y,t_min,t_max,i1,i2,ave_y,err_y
  
     ave_str   = strmid(strtrim(string(ave_y),2),0,5)
     err_str   = strmid(strtrim(string(err_y),2),0,5)

  ;;    trend_str = strmid(strtrim(string(res.trend),2),0,5)
  ;;    hump_str  = strmid(strtrim(string(res.hump),2),0,5)

     text = ' !3['+ave_str+'!9+!3'+err_str+']'
   
     ;; screen prints  ave_y, ave_y1, ave_y2, ave_y3

      print, 'ave_y=',ave_y
      print, 'err_y=',err_y
 

  ;;         res = meanstddev(pname,0,y1,t,t_min,t_max)
  ;;         ave_y1 = res.avg
  ;;         err_y1 = res.sig2

       diff_statistics,y1,t_min,t_max,i1,i2,ave_y1,err_y1

      print, 'ave_y1=',ave_y1
      print, 'err_y1=',err_y1

  ;;          res = meanstddev(pname,0,y2,t,t_min,t_max)
  ;;         ave_y2 = res.avg
  ;;         err_y2 = res.sig2

       diff_statistics,y2,t_min,t_max,i1,i2,ave_y2,err_y2

      print, 'ave_y2=',ave_y2
      print, 'err_y2=',err_y2

 ;;        res = meanstddev(pname,0,y3,t,t_min,t_max)
 ;;        ave_y3 = res.avg
 ;;        err_y3 = res.sig2

       diff_statistics,y3,t_min,t_max,i1,i2,ave_y3,err_y3

      print, 'ave_y3=',ave_y3
      print, 'err_y3=',err_y3

     ;; end screen prints


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
     
     ;; Diffusion trace

     oplot,t,y,color=color_vec[0]
     ;; 0=blue


     ;; Average line

     oplot,t[it1:it2],unity*ave_y,color=color_vec[1]   

     ;; Curve dots

     oplot,t[it1]*[1,1],y[it1]*[1,1],psym=8,color=line   
     oplot,t[it2]*[1,1],y[it2]*[1,1],psym=8,color=line   

     ;; RMS deviation bars

     oplot,t[it1:it2],unity*(ave_y-err_y),color=color_vec[2]   
     oplot,t[it1:it2],unity*(ave_y+err_y),color=color_vec[2]   

    if i_tp eq 0 then begin

      oplot,t,y1,color=color_vec[1]
  ;; 1=purple  M=z H=par
      oplot,t,y2,color=color_vec[2]
  ;; 2=yellow  M=y H=perp
      oplot,t,y3,color=color_vec[3]
  ;; 3=green   M=u H=check

     oplot,t,y4,color=color_vec[0]

     endif


     ynorm = ymax/zoom

   ;;  xyouts,0.7*max(t),0.95*ynorm,$
   ;;    '!3trend='+trend_str

   ;;  xyouts,0.7*max(t),0.9*ynorm,$
   ;;    '!3hump='+hump_str

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
     
  endif else begin

     plot_def_new,pname+'_PDF'

     pdf_statistics,y,it1,it2,y_bin,pdf

     y_bin_plot = fltarr(2*nbin)
     pdf_plot   = fltarr(2*nbin)

     i  = indgen(nbin)
     dy = y_bin[1]-y_bin[0]

     y_bin_plot[2*i]   = y_bin[i]-0.5*dy
     y_bin_plot[2*i+1] = y_bin[i]+0.5*dy
     pdf_plot[2*i]     = pdf[i]
     pdf_plot[2*i+1]   = pdf[i]

     xmin = min(y_bin_plot)
     xmax = max(y_bin_plot)
     ymin = min(pdf_plot)
     ymax = max(pdf_plot)

     yave = total(y_bin[*]*pdf[*])

     plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[xmin,xmax],$
       xtitle=ytitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax],$
       ytitle='!3Probability Density Function', $
       color=line

     oplot,y_bin_plot,pdf_plot,color=color_vec[0]
     oplot,yave*[1,1],100*[-1,1],color=color_vec[0],linestyle=1

  endelse

  plot_finish

  return

end
