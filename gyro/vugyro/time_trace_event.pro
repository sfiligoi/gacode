pro time_trace_event, time_trace

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_TIME_TRACE
  
  widget_control, time_trace.id, $
    get_uvalue=uvalue

  wset, widget
  
  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: begin
        index=1
        goto, plot_it
     end
     2: begin
        index=2
        goto, plot_it
     end
     3: begin
        index=3
        goto, plot_it
     end
     4: begin
        index=4
        goto, plot_it
     end

     5: begin
        counter_up,j_c,n_theta_plot-1,1
        goto, plot_it
     end

     6: begin
        counter_dn,j_c,0,1 
        goto, plot_it
     end

     7: begin
        counter_up,i_c,n_r-1,1
        goto, plot_it
     end

     8: begin
        counter_dn,i_c,0,1 
        goto, plot_it
     end

     81: begin
        i_er_shear=i_er_shear+1
        if(i_er_shear eq 3) then begin
           i_er_shear = 0
        endif
        goto, plot_it
     end

     9: begin
        counter_up,in_c,n_n-1,1
        goto, plot_it
     end

     10: begin
        counter_dn,in_c,0,1 
        goto, plot_it
     end

     11: begin
        i_log = 1-i_log 
        goto, plot_it
     end

     12: begin
        counter_up,t_dot,n_time1,1 
        goto, plot_it
     end

     13: begin
        counter_dn,t_dot,0,1
        goto, plot_it
     end

     14: begin
        if (n_field gt 1) then i_f = 1-i_f
        goto, plot_it
     end

     15: begin
        plot_mode = 2
        goto, plot_it
     end
     
     16: begin
        plot_export = 1
        goto, plot_it
     end

     17: widget_control, time_trace.top, /destroy

  endcase
  
  ;;----------------------------------------------------------

  return

  plot_it:

  time_trace_plot

end
