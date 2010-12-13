pro spectrum_np_event, spectrum_np

  common GLOBAL
  common PRIVATE_SPECTRUM_NP
  
  widget_control, spectrum_np.id, $
                  get_uvalue=uvalue

  wset, widget
  
  ;;-------------------------------------------------------
  ;; Now, begin menu selection
  ;;
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
        counter_up,t_c,n_time1,20 
        goto, plot_it
     end

     5: begin 
        counter_dn,t_c,0,20
        goto, plot_it
     end

     6: begin
        counter_up,ax_c,90,5
        goto, plot_it
     end

     7: begin
        p_squ = 1-p_squ
        goto, plot_it
     end

     8: begin
        counter_dn,ax_c,-90,5
        goto, plot_it
     end

     9: begin
        counter_up,az_c,90,5
        goto, plot_it
     end

     10: begin 
        counter_dn,az_c,-90,5
        goto, plot_it
     end

     11: begin 
        i_log = 1-i_log
        goto, plot_it
     end

     12: begin
        plot_mode = 2
        goto, plot_it
     end

     13: widget_control, spectrum_np.top, /destroy

  endcase                              

  return

  plot_it:

  widget_control, spectrum_np.top, $
                  get_uvalue=state, $
                  /no_copy
  
  t_str = 't = '+strtrim(string(t[t_c]),2)

  widget_control, state.t_label, set_value=t_str

  widget_control, spectrum_np.top, $
                  set_uvalue=state, $
                  /no_copy
  
  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  spectrum_np_plot
  
  return

end
