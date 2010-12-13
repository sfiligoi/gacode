pro spectrum_np_ave_event, spectrum_np_ave

  common GLOBAL
  common PRIVATE_SPECTRUM_NP_AVE
  
  widget_control, spectrum_np_ave.id, $
                  get_uvalue=uvalue

  wset, widget
  
  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        counter_up,ax_c,90,5
        goto, plot_it
     end

     3: begin
        counter_dn,ax_c,-90,5
        goto, plot_it
     end

     4: begin
        counter_up,az_c,90,5
        goto, plot_it
     end

     5: begin 
        counter_dn,az_c,-90,5
        goto, plot_it
     end

     6:  begin
        p_squ = 1-p_squ
        goto, plot_it
     end

     7: begin 
        i_log = 1-i_log
        goto, plot_it
     end

     8: begin
        plot_mode = 2
        goto, plot_it
     end

     9: widget_control, spectrum_np_ave.top, /destroy
     
  endcase                              
  
  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------
  
  spectrum_np_ave_plot

  return

end
