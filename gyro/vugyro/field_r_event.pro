pro field_r_event, field_r

  common GLOBAL

  widget_control, field_r.id, $
                  get_uvalue=uvalue

  wset, field_r_wid

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        i_all = 1-i_all
        goto, plot_it
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

     11: begin
        if (n_field gt 1) then i_f = 1-i_f
        goto, plot_it
     end

     12: begin
        i_abs = 1-i_abs
        goto, plot_it
     end

     13: begin
        i_log = 1-i_log
        goto, plot_it
     end

     14: begin
        plot_mode = 2
        goto, plot_it
     end

     15: begin
        plot_export = 1
        goto, plot_it
     end

     16: widget_control, field_r.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; WIDGET LABELS
  ;;-------------------------------------------------------

  widget_control, field_r.top, $
                  get_uvalue=state, $
                  /no_copy

  j_str = ' theta/pi = '+strtrim(string(thetad_plot[j_c]),2)
  
  widget_control, state.j_label, set_value=j_str

  widget_control, field_r.top, $
                  set_uvalue=state, $
                  /no_copy

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  field_r_plot

  return

end
