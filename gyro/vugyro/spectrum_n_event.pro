pro spectrum_n_event, spectrum_n

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_SPECTRUM_N

; In auto mode  spectrum_n won't be defined
if(n_elements(spectrum_n) le 0) then goto,plot_it

  widget_control, spectrum_n.id, $
    get_uvalue=uvalue

  wset, widget

;;-------------------------------------------------------
;; MENU
;;-------------------------------------------------------

  case (uvalue) of 
     
     1: begin
        active_spectrum_n = 1  
        title  = '!3RMS n = '
        ytitle = '!4<u>!3(t)'
        pname0  = 'spectrum_n_'
        goto, plot_it
     end   

     2: begin
        active_spectrum_n = 2
        title  = '!3RMS Zonal and finite-n potentials'
        ytitle = '!4<u>!3(t)'
        pname0  = 'spectrum_n_all'
        goto, plot_it
     end   

     3: begin
        counter_up,in_c,n_n-1,1
        goto, plot_it
     end

     4: begin
        counter_dn,in_c,0,1
        goto, plot_it
     end

     5: begin
        zoom = zoom*2
        goto, plot_it
     end

     6: begin
        zoom = zoom/2
        goto, plot_it
     end

     7: begin
        plot_mode = 2
        goto, plot_it
     end

     8: begin
        plot_export = 1
        goto, plot_it
     end

     9: widget_control, spectrum_n.top, /destroy

  endcase                              

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  case (active_spectrum_n) of

     0: return

     1: begin

        pname = pname0+strtrim(string(n_tor[in_c]),2)
        plot_def_new,pname

        y = fltarr(n_time)
        for tt=0,n_time1 do begin
           ;; Removed 1/n_r because a sum over p gives the
           ;; proper radial average of phi^2.
           y[tt] = total(kxkyspec[*,in_c,tt])
        endfor

        y = sqrt(y)
        y_axis_magic,y,ymin,ymax,d_y

        plot,[0],[1],$
          /nodata,$
          title=title+strtrim(string(n_tor[in_c]),2)+' potential',$
          xstyle=1,$
          xminor=0,$
          xrange=[min(t),max(t)],$
          xtitle=csa_string,$
          ystyle=1,$
          yminor=0,$
          yrange=[0,ymax]/zoom,$
          color=line

        oplot,t,y,color=color_vec[0]   

        ;;---------------------------------------------------
        ;; DATA EXPORT
        ;;
        if (plot_export eq 1) then begin
           openw,1,pname+'.idlout'
           for tt=0,n_time-1 do begin
              printf,1,t[tt],y[tt]
           endfor
           close,1
           print,'Exported data to ',pname+'.idlout'
           plot_export = 0
        endif
        ;;---------------------------------------------------

     end

     2: begin

        if (n_n eq 1) then return

        pname = pname0
        plot_def_new,pname

        y = fltarr(n_n,n_time)
        for nn=0,n_n-1 do begin
           for tt=0,n_time1 do begin
              y[nn,tt] = total(kxkyspec[*,nn,tt])
           endfor
        endfor

        yn = fltarr(2,n_time)
        
        for tt=0,n_time1 do begin
           yn[1,tt] = sqrt(2.0*total(y[1:n_n-1,tt]))
           yn[0,tt] = sqrt(y[0,tt])
        endfor

;  y_axis_magic is for arrays with a single dimension, so this doesn't work:
        y_axis_magic,yn,ymin,ymax,d_y
        ymax=max(yn) ; this looks at all of yn

        spectrum_n_max = ymax/zoom

        plot,[0],[1],$
          /nodata,$
          title=title,$
          ytitle=ytitle,$
          xstyle=1,$
          xminor=0,$
          xrange=[0.0,max(t)],$
          xtitle=csa_string,$
          ystyle=1,$
          yminor=0,$
          yrange=spectrum_n_max*[0,1],$
          color=line

        oplot,t,yn[0,*],color=color_vec[0],linestyle=line_vec[0] 
        oplot,t,yn[1,*],color=color_vec[1],linestyle=line_vec[1]    

        xyouts,0.1*max(t),$
          0.95*spectrum_n_max,'!3n=0',$
          color=color_vec[0]

        xyouts,0.1*max(t),$
          0.9*spectrum_n_max,'!3n>0',$
          color=color_vec[1]

        ;;---------------------------------------------------
        ;; DATA EXPORT
        ;;
        if (plot_export eq 1) then begin
           openw,1,pname+'.idlout'
           for tt=0,n_time-1 do begin
              printf,1,t[tt],yn[0,tt],yn[1,tt]
           endfor
           close,1
           print,'Exported data to ',pname+'.idlout'
           plot_export = 0
        endif
        ;;---------------------------------------------------

     end

  endcase

  plot_finish

  return

end
