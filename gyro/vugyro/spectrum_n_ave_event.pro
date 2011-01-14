pro spectrum_n_ave_event, spectrum_n_ave

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_SPECTRUM_N_ave

  widget_control, spectrum_n_ave.id, $
                  get_uvalue=uvalue

  wset, widget

;;-------------------------------------------------------
;; MENU
;;-------------------------------------------------------

  case (uvalue) of 
     
     1: begin
        goto, plot_it
     end   

     2: begin
        zoom = zoom*2
        goto, plot_it
     end

     3: begin
        zoom = zoom/2
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

     6: widget_control, spectrum_n_ave.top, /destroy

  endcase                              

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  title  = 'Radial asymmetry'
  ytitle = '<'+kr_rho_string+'!4u!u2!n>/<u!u2!n>!3'
  pname  = 'spectrum_n_ave_'
  plot_def_new,pname

  yt1 = fltarr(n_time)
  yt2 = fltarr(n_time)
  y  = fltarr(n_n)
  kxrho = (indgen(n_r)-n_r/2)*2*!pi/(r[n_r-1]-r[0])*rho_s
  for in=0,n_n-1 do begin
      yt1[*] = 0.0
      yt2[*] = 0.0
      for ir=0,n_r-1 do begin
         yt1[*] = yt1[*]+kxrho[ir]*kxkyspec[ir,in,*]
         yt2[*] = yt2[*]+kxkyspec[ir,in,*]
      endfor
      diff_stat_fast,yt1,it1,it2,ave_y1
      diff_stat_fast,yt2,it1,it2,ave_y2
      y[in] = ave_y1/ave_y2
  endfor

  plot,[0],[1],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[0.0,max(kt_rho)],$
       xtitle=kt_rho_string,$
       ystyle=1,$
       yminor=0,$
       yrange=max(abs(y))*[-1,1]/zoom,$
       ytitle=ytitle,$
       color=line

  oplot,kt_rho,y,color=color_vec[0]   

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for tt=0,n_time-1 do begin
        printf,1,kr_rho[tt],y[tt]
     endfor
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
