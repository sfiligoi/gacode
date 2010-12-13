;; plot radial shear in (phi,A)

pro field_shear_event, field_shear

  common GLOBAL
  common PRIVATE_FIELD_SHEAR
  
  widget_control, field_shear.id, $
    get_uvalue=uvalue

  wset, widget

  xtitle = '!3r/a'
  if i_f eq 0 then begin
     title  = '!3d!u2!n!4u/!3dr!u2!n'
  endif else begin 
     title  = '!3d!u2!n!4A/!3dr!u2!n'
  endelse


  ;;-------------------------------------------------------
  ;; Now, begin menu selection
  ;;
  top:

  case (uvalue) of 
     
     1: begin

        active_f = 1

        ;;-------------------------------------------------
        ;; widget control for label functions
        ;;
        widget_control, field_shear.top, $
          get_uvalue=state, $
          /no_copy
        
        j_str = 'theta/pi = '+strtrim(string(thetad_plot[j_c]),2)
        n_str = 'n = '+strtrim(string(n_tor[in_c]),2)

        widget_control, state.t_label, set_value=t_string
        widget_control, state.j_label, set_value=j_str
        widget_control, state.n_label, set_value=n_str      

        widget_control, field_shear.top, $
          set_uvalue=state, $
          /no_copy
        ;;-------------------------------------------------

        y0 = fltarr(2,n_r,n_time)
        yp = fltarr(2,n_r,n_time)

        y0[*,*,*] = rho_s*U[*,j_c,*,i_f,in_c,*]
        for i=2,n_r-3 do begin
           yp[*,i,*] = (-(1.0/12.0)*y0[*,i+2,*] $
                        +(16.0/12.0)*y0[*,i+1,*] $
                        -(30.0/12.0)*y0[*,i,*] $ 
                        +(16.0/12.0)*y0[*,i-1,*] $ 
                        -(1.0/12.0)*y0[*,i-2,*])/ $
             (r[1]-r[0])^2
        endfor

        yp[*,0,*]     = yp[*,2,*]
        yp[*,1,*]     = yp[*,2,*]
        yp[*,n_r-1,*] = yp[*,n_r-3,*]
        yp[*,n_r-2,*] = yp[*,n_r-3,*]

        ymax = 1.1*max(abs(yp)) 

        plot_def_new,'field_shear_all'

        plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[min(r),max(r)],$
          xtitle=xtitle,$
          ystyle=1,$
          yminor=0,$
          yrange=[-1,1]*ymax,$
          ytitle='!3aqua=Re ; yellow=Im',$
          color=line

        for tt=0,n_time1 do begin
           if (t[tt] gt t_min) and (t[tt] lt t_max) then begin
              oplot,r,yp[0,*,tt],$
                color=color_vec[0],$
                linestyle=line_vec[0]
           endif
        endfor

        plot_finish

     end   

     2: begin

        active_f = 2

        ;;-------------------------------------------------
        ;; widget control for label functions
        ;;
        widget_control, field_shear.top, $
          get_uvalue=state, $
          /no_copy
        
        t_str = 't = '+strtrim(string(t[t_c]),2)
        j_str = 'theta/pi = '+strtrim(string(thetad_plot[j_c]),2)
        n_str = 'n = '+strtrim(string(n_tor[in_c]),2)

        widget_control, state.t_label, set_value=t_str
        widget_control, state.j_label, set_value=j_str
        widget_control, state.n_label, set_value=n_str      

        widget_control, field_shear.top, $
          set_uvalue=state, $
          /no_copy
        ;;-------------------------------------------------


        y0 = fltarr(2,n_r,n_time)
        yp = fltarr(2,n_r,n_time)

        y0[*,*,*] = rho_s*U[*,j_c,*,i_f,in_c,*]
        for i=2,n_r-3 do begin
           yp[*,i,*] = (-(1.0/12.0)*y0[*,i+2,*] $
                        +(16.0/12.0)*y0[*,i+1,*] $
                        -(30.0/12.0)*y0[*,i,*] $ 
                        +(16.0/12.0)*y0[*,i-1,*] $ 
                        -(1.0/12.0)*y0[*,i-2,*])/ $
             (r[1]-r[0])^2
        endfor

        yp[*,0,*]     = yp[*,2,*]
        yp[*,1,*]     = yp[*,2,*]
        yp[*,n_r-1,*] = yp[*,n_r-3,*]
        yp[*,n_r-2,*] = yp[*,n_r-3,*]

        ymax = 1.1*max(abs(yp)) 

        plot_def_new,'field_shear'

        plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[min(r),max(r)],$
          xtitle=xtitle,$
          ystyle=1,$
          yminor=0,$
          yrange=[-1,1]*ymax,$
          ytitle='!3aqua=Re ; yellow=Im',$
          color=line

        oplot,r,yp[0,*,t_c],$
          color=color_vec[0],$
          linestyle=line_vec[0]    
        
        oplot,r,yp[0,*,t_c],color=line,psym=8,symsize=0.5

        oplot,r,0.0*r,color=line,linestyle=1

        plot_finish

     end   
     
     3: begin
        counter_up,t_c,n_time1,1
        goto, plot_it
     end

     4: begin
        counter_dn,t_c,0,1
        goto, plot_it
     end

     5: begin
        counter_up,t_c,n_time1,10 
        goto, plot_it
     end

     6: begin 
        counter_dn,t_c,0,10
        goto, plot_it
     end

     7: begin
        counter_up,in_c,n_n-1,1 
        goto, plot_it
     end

     8: begin
        counter_dn,in_c,0,1  
        goto, plot_it
     end

     9: begin
        counter_up,j_c,n_theta_plot-1,1  
        goto, plot_it
     end

     10: begin
        counter_dn,j_c,0,1
        goto, plot_it
     end

     13: begin
        if (n_field gt 1) then i_f = 1-i_f
        goto, plot_it
     end

     14: widget_control, field_shear.top, /destroy

     15: begin
        plot_mode = 2
        uvalue = active_f
        goto,top
     end

  endcase
  ;;
  ;;----------------------------------------------------------

  return

  plot_it:
  
  uvalue=active_f
  goto,top

end
