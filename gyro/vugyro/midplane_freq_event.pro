pro midplane_freq_event, midplane_freq

  common GLOBAL
  common MIDPLANE_DATA
  common PLOT_VARIABLES
  common PRIVATE_MIDPLANE_FREQ

  widget_control, midplane_freq.id, $
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
        counter_up,i_pwr,n_pwr-1,1
        goto, plot_it
     end   

     3: begin
        counter_dn,i_pwr,0,1
        goto, plot_it
     end   

     5: begin
        log_flag = 1-log_flag
        goto, plot_it
     end

     6: begin
        spec_amp = spec_amp*2.0
        goto, plot_it
     end

     7: begin
        spec_amp = spec_amp/2.0
        goto, plot_it
     end

     8: begin
        i_ptype = 1-i_ptype
        goto, plot_it
     end

     9: begin
        plot_mode = 2
        goto, plot_it
     end

     10: begin
        plot_export = 1
        goto, plot_it
     end

     11: widget_control, midplane_freq.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------
  
  nt = it2-it1+1

  if (log_flag eq 0) then begin
     logstr = ' (linear scale)'
  endif else begin
     logstr = ' (log scale)' 
  endelse

  title = tag[i_pwr]+' frequency response'
  pname = 'tspec_'+tag[i_pwr]

  yt = complexarr(nt)
  z = fltarr(nt,n_n)

  w = -2*!pi/(t_max-t_min)*(indgen(nt)-(nt/2-1)) 

  xtitle = '!3(a/c!ds!n)!4x!3 for'+t_string
  ytitle = 'Toroidal mode number'+logstr

  for i_n=0,n_n-1 do begin

     ;; Hardwired central radius -> n_r/2
     
     yt[*] = pwr[i_pwr,n_r/2,i_n,it1:it2]
     z[*,i_n] = shift(abs(fft(yt[*])),nt/2-1)

  endfor

  y_ax = fltarr(n_n+1)
  for i_n=0,n_n do begin
     y_ax[i_n] = n_dn*(i_n-0.5)
  endfor

  if (log_flag eq 1) then z = alog10(z)
  
  if (i_ptype eq 0) then begin

     plot_def_new,pname+'_CMAP'

     z_min = min(z)
     z_max = max(z)

     clevel = z_min+findgen(nlevels)*(z_max-z_min)/(nlevels-1)

     loadct,CTab  

     for i_n=0,n_n-1 do begin

        yt[*] = z[*,i_n]

        if (i_n eq 0) then begin

           contour,yt#[1,1],w,[y_ax[i_n],y_ax[i_n+1]],$
             ystyle=1,yrange=[min(y_ax),max(y_ax)],$
             xstyle=1,xrange=[-4.0,4.0],$
             xtitle=xtitle,$
             ytitle=ytitle,$
             title='!3'+title,$
             nlevels=nlevels,$
             levels=clevel,$
             /fill

        endif else begin

           contour,yt#[1,1],w,[y_ax[i_n],y_ax[i_n+1]],$
             ystyle=1,yrange=[min(y_ax),max(y_ax)],$
             xstyle=1,xrange=[-1.5,1.5],$
             xtitle=xtitle,$
             ytitle=ytitle,$
             title='!3'+title,$
             nlevels=nlevels,$
             levels=clevel,$
             /fill,$
             /noerase

        endelse

     endfor

     set_line_colors

     for i_n=0,n_n-1 do begin
        oplot,10*[-1,1],y_ax[i_n+1]*[1,1]
     endfor

     oplot,[0,0],[min(y_ax),max(y_ax)],linestyle=1,color=color_vec[1]

  endif else begin

     plot_def_new,pname+'_LINE'

     z_max = max(z)

     plot,[0],[0],$
       /nodata, $
       ytitle=ytitle, $
       ystyle=1,yrange=[min(y_ax),max(y_ax)],$
       xtitle=xtitle, $
       xstyle=1,xrange=[-1.5,1.5],$
       title='!3'+title

     for i_n=0,n_n-1 do begin
        oplot,w,spec_amp*2*z[*,i_n]/z_max+i_n*n_dn     
     endfor

     oplot,[0,0],[min(y_ax),max(y_ax)],linestyle=1,color=color_vec[1]

  endelse

  plot_finish

  return

end
