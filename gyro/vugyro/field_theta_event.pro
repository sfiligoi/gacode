;; plot theta-dependence of (phi,A)

pro field_theta_event, field_theta

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_FIELD_THETA
  
  widget_control, field_theta.id, $
      get_uvalue=uvalue

  wset, widget

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
      widget_control, field_theta.top, $
          get_uvalue=state, $
          /no_copy
      
      t_str = 't = 0-'+strtrim(string(t[n_time1]),2)
      r_str = 'r = '+strtrim(string(r[i_c]),2)
      n_str = 'n = '+strtrim(string(n_tor[in_c]),2)

      widget_control, state.t_label, set_value=t_str
      widget_control, state.r_label, set_value=r_str
      widget_control, state.n_label, set_value=n_str      

      widget_control, field_theta.top, $
          set_uvalue=state, $
          /no_copy
      ;;-------------------------------------------------

      ymax = 1.1*max(abs(u[0,*,i_c,i_f,in_c,*])) 

      if i_f eq 0 then begin
        plot_def_new,'phi_theta_'+strtrim(string(fix(n_tor[in_c])),2)
        title = '!4u(h)'
      endif else begin 
        plot_def_new,'a_theta_'+strtrim(string(fix(n_tor[in_c])),2)
        title = '!3A!d!9#!4!n(h)'
      endelse

      plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[-1,1],$
          xtitle='!4h!3!dj!n/!4p!3',$
          ystyle=1,$
          yminor=0,$
          yrange=[-1,1]*ymax,$
          ytitle='!3aqua=Re ; yellow=Im',$
          color=line

      for tt=0,n_time1 do begin
        oplot,thetad_plot,u[0,*,i_c,i_f,in_c,tt],color=color_vec[0]
      endfor
      for tt=0,n_time1 do begin
        oplot,thetad_plot,u[1,*,i_c,i_f,in_c,tt],color=color_vec[1]
      endfor

      plot_finish

    end   

    2: begin

      active_f = 2

      ;;-------------------------------------------------
      ;; widget control for label functions
      ;;
      widget_control, field_theta.top, $
          get_uvalue=state, $
          /no_copy
      
      t_str = 't = '+strtrim(string(t[t_c]),2)
      r_str = 'r = '+strtrim(string(r[i_c]),2)
      n_str = 'n = '+strtrim(string(n_tor[in_c]),2)

      widget_control, state.t_label, set_value=t_str
      widget_control, state.r_label, set_value=r_str
      widget_control, state.n_label, set_value=n_str      

      widget_control, field_theta.top, $
          set_uvalue=state, $
          /no_copy
      ;;-------------------------------------------------

      ymax = 1.1*max(abs(u[*,*,i_c,i_f,in_c,t_c])) 

      if i_f eq 0 then begin
        plot_def_new,'phi_theta_t_'+strtrim(string(fix(n_tor[in_c])),2)
        title = '!4u(h)'
      endif else begin 
        plot_def_new,'a_theta_t_'+strtrim(string(fix(n_tor[in_c])),2)
        title = '!3A!d!9#!4!n(h)'
      endelse

      plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[-1,1],$
          xtitle='!4h!3!dj!n/!4p!3',$
          ystyle=1,$
          yminor=0,$
          yrange=[-1,1]*ymax,$
          ytitle='!3aqua=Re ; yellow=Im',$
          color=line

      oplot,thetad_plot,u[0,*,i_c,i_f,in_c,t_c],color=color_vec[0]    
      oplot,thetad_plot,u[1,*,i_c,i_f,in_c,t_c],color=color_vec[1]  
      oplot,thetad_plot,u[0,*,i_c,i_f,in_c,t_c],color=line,psym=8   
      oplot,thetad_plot,u[1,*,i_c,i_f,in_c,t_c],color=line,psym=8

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
      counter_up,i_c,n_r-1,1 
      goto, plot_it
    end

    10: begin
      counter_dn,i_c,0,1
      goto, plot_it
    end

    11: begin
      if (n_field gt 1) then i_f = 1-i_f
      goto, plot_it
    end

    12: widget_control, field_theta.top, /destroy
    
    13: begin
      plot_mode = 2
      uvalue = active_f
      goto,top
    end

  endcase
  
  return

  plot_it:
  
  uvalue=active_f
  goto,top

end
