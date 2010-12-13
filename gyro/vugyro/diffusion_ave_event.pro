pro diffusion_ave_event, diffusion_ave

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_DIFFUSION_AVE  

  widget_control, diffusion_ave.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     0: goto, plot_it

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
        i_f = i_f+1
        if (i_f ge n_field) then i_f = -1
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
        i_units = 1-i_units
        goto, plot_it
     end

     7: begin
        if (exists_diff_t eq 1) then begin
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

     10: widget_control, diffusion_ave.top, /destroy

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
  
  smf_tag,title,pname,ytitle

  case (i_tp) of

     0: begin

        pname = 'diff-'+pname
        title = title+' Diffusion'
        if i_f ge 0 then begin
           y[*] = diff[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = diff[i_spec,0,i_moment,*]+$
             diff[i_spec,1,i_moment,*]
        endelse

     end

     1: begin

        pname = 'diff_t-'+pname
        title = title+' (trap) Diffusion'
        if i_f ge 0 then begin
           y[*] = diff_trapped[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = diff_trapped[i_spec,0,i_moment,*]+$
             diff_trapped[i_spec,1,i_moment,*]
        endelse

     end

     2: begin

        pname = 'diff_p-'+pname
        title = title+' (pass) Diffusion'
        if i_f ge 0 then begin
           y[*] = diff[i_spec,i_f,i_moment,*]-$
             diff_trapped[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = diff[i_spec,0,i_moment,*]-$
             diff_trapped[i_spec,0,i_moment,*]+$
             diff[i_spec,1,i_moment,*]-$
             diff_trapped[i_spec,1,i_moment,*]
        endelse

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

     y_axis_magic,y,ymin,ymax,d_y

     unity = fltarr(it2-it1)
     unity[*] = 1.0

     if (it2-it1) gt 10 then begin

        ;; Enough points for Mikkelsen's method

        res = meanstddev(pname,0,y,t,t_min,t_max)
        ave_y = res.avg
        err_y = res.sig2
        trend_y = res.trend
        hump_y = res.hump

     endif else begin

        ;; Only compute simple average

        diff_stat_fast,y,it1,it2,ave_y
        err_y = -1
        trend_y = -1
        hump_y = -1

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
     
     ;; Diffusion trace

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
     
  endif else begin

     ;; Diffusion PDF histogram

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
