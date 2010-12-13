pro linear_freq_event, linear_freq

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_LINEAR_FREQ
  
  widget_control, linear_freq.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it
     
     2: begin
        counter_up,in_c,n_n-1,1
        goto, plot_it
     end

     3: begin
        counter_dn,in_c,0,1 
        goto, plot_it
     end

     4: begin
        real_flag = 1-real_flag
        goto, plot_it
     end

     5: begin
        zoom = zoom/sqrt(2)
        goto,plot_it
     end

     6: begin
        zoom = zoom*sqrt(2)
        goto,plot_it
     end

     7: begin
        plot_mode = 2
        goto, plot_it
     end

     8: widget_control, linear_freq.top, /destroy
     
  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  title = kt_rho_string+'='+strmid(strtrim(string(kt_rho[in_c]),2),0,5)

  if (real_flag eq 0) then begin
     ytitle = '!3(a/c!ds!n) !4x!3(t)'
     pname = 'linfreq_re'
  endif else begin
     ytitle = '!3(a/c!ds!n) !4c!3(t)'
     pname = 'linfreq_im'
  endelse

  plot_def_new,pname
  
  wn = fltarr(n_time)

  wn[*] = w[real_flag,in_c,*]
  diff_stat_fast,wn,it1,it2,ave_w
  
  ave_str = strmid(strtrim(string(ave_w),2),0,6)

  plot,[0],[0],$
    /nodata,$
    title=title+' ['+ave_str+']',$
    xstyle=1,$
    xrange=[0,max(t)],$
    xtitle=csa_string,$
    ystyle=1,$
    yrange=[ave_w-zoom,ave_w+zoom],$
    ytitle=ytitle,$
    color=line

  ;; time trace  
  oplot,t,w[real_flag,in_c,*],color=color_vec[0]
  
  ;; average over [t_min,t_max]
  oplot,[0,max(t)],ave_w*[1,1],linestyle=1
  oplot,[t_min,t_max],ave_w*[1,1]

  plot_finish

end
