pro linear_spec_event, linear_spec_obj

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_LINEAR_SPEC
  

  widget_control, linear_spec_obj.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it
     
     2: begin
        real_flag = 1-real_flag
        goto, plot_it
     end

     3: begin
        zero_flag = 1-zero_flag
        goto, plot_it
     end

     4: begin
        plot_mode = 2
        goto, plot_it
     end

     5: begin
        plot_export = 1
        goto, plot_it
     end

     6: widget_control, linear_spec_obj.top, /destroy
     
  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  if (real_flag eq 0) then begin
     ytitle = '!3(a/c!ds!n) !4x'
     pname = 'linspec_re'
  endif else begin
     ytitle = '!3(a/c!ds!n) !4c'
     pname = 'linspec_im'
  endelse

  xtitle = kt_rho_string+'!3 with'+t_string

  plot_def_new,pname
  
  w0 = fltarr(n_n,n_time)
  w0[*,*] = w[real_flag,*,*]
  wt = fltarr(n_time)
  wn = fltarr(n_n)

  for i_n=0,n_n-1 do begin
     wt[*] = w0[i_n,*]  
     diff_stat_fast,wt,it1,it2,ave_w
     wn[i_n] = ave_w
  endfor

  ymin = 1.1*min(wn)
  if (ymin gt 0.0) then ymin = 0.0

  ymax = 1.1*max(wn)
  if (ymax lt 0.0) then ymax = 0.0

  if (n_n gt 1) then begin
     dk = kt_rho[1]-kt_rho[0]
     xmin = (kt_rho[0]-0.5*dk)*(1-zero_flag)
     xmax = kt_rho[n_n-1]+0.5*dk
  endif else begin
     xmin = 0.0
     xmax = 2*kt_rho[0]
  endelse

  plot,[1],[1],$
    /nodata,$
    title=title,$
    xstyle=1,$
    xminor=0,$
    xrange=[xmin,xmax],$
    xtitle=xtitle,$
    ystyle=1,$
    yminor=0,$
    yrange=[ymin,ymax],$
    ytitle=ytitle,$
    color=line

  oplot,kt_rho,wn,color=color_vec[0]
  oplot,kt_rho,wn,color=line,psym=8

  ;; Attach to the origin if flag set:
  if zero_flag eq 1 then begin
     oplot,[0,kt_rho[0]],[0,wn[0]],color=line
  endif

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i_n=1,n_n-1 do begin
        printf,1,kt_rho[i_n],wn[i_n]
     endfor
     printf,1,' Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
