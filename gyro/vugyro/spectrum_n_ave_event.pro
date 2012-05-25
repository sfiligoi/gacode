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
  pname  = 'spectrum_n_ave'
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
     for in=0,n_n-1 do begin
        printf,1,kt_rho[in],y[in]
     endfor
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  ;; output time average kxky spectrum to file
     t_total = t[it2]-t[it1]
     sumkxky=fltarr(n_r,n_n)
     avekxky=fltarr(n_r,n_n)
     sumkxky[*,*] = 0.0
     for i=it1,it2-1 do begin
       dt = t[i+1]-t[i]
       sumkxky[*,*] = sumkxky[*,*]+0.5*(kxkyspec[*,*,i+1]+kxkyspec[*,*,i])*dt
    endfor
    avekxky[*,*]=sumkxky[*,*]/t_total
    openw,1,'avekxky.idlout'
    kxrho = (indgen(n_r)-n_r/2)*2*!pi/(r[n_r-1]-r[0])*rho_s
    printf,1,'kx',kt_rho[0:n_n-1],format='(a3,128(2x,e15.5))'
    for ir=0,n_r-1 do begin
        printf,1,kxrho[ir],avekxky[ir,0:n_n-1],format='(128(2x,e15.5))'
    endfor
    close,1
    print,'Exported time average kxky spectrum to avekxky.idlout'
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
