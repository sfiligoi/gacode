;; plot theta-dependence of (phi,A) at r=r0

pro field_r0_theta_event, field_r0_theta

  common GLOBAL
  common PLOT_VARIABLES
  
  widget_control, field_r0_theta.id, $
    get_uvalue=uvalue

  wset, field_r0_wid

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     1: goto, plot_it
     
     2: begin
        counter_up,t_c,n_time1,1
        goto, plot_it
     end

     3: begin
        counter_dn,t_c,0,1
        goto, plot_it
     end

     4: begin
        counter_up,t_c,n_time1,10 
        goto, plot_it
     end

     5: begin 
        counter_dn,t_c,0,10
        goto, plot_it
     end

     6: begin
        counter_up,in_c,n_n-1,1 
        goto, plot_it
     end

     7: begin
        counter_dn,in_c,0,1  
        goto, plot_it
     end

     8: begin
        if (n_field gt 1) then i_f = 1-i_f
        goto, plot_it
     end

     
     9: begin
        plot_mode = 2
        goto, plot_it
     end

     10: widget_control, field_r0_theta.top, /destroy

  endcase
  
  return

  plot_it:
  
  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  ymax = 1.1*max(abs(field_r0[*,*,i_f,in_c,t_c])) 

  if i_f eq 0 then begin
     plot_def_new,'phi_r0_t_'+strtrim(string(fix(n_tor[in_c])),2)
     title = '!4u(h)'
  endif else begin 
     plot_def_new,'a_r0_t_'+strtrim(string(fix(n_tor[in_c])),2)
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

  oplot,thetad_r0_plot,field_r0[0,*,i_f,in_c,t_c],color=color_vec[0]    
  oplot,thetad_r0_plot,field_r0[1,*,i_f,in_c,t_c],color=color_vec[1]  
  oplot,thetad_r0_plot,field_r0[0,*,i_f,in_c,t_c],color=line,psym=8   
  oplot,thetad_r0_plot,field_r0[1,*,i_f,in_c,t_c],color=line,psym=8

  plot_finish

  return

end
