pro midplane_power_event, midplane_power

  common GLOBAL
  common MIDPLANE_DATA

  widget_control, midplane_power.id, $
    get_uvalue=uvalue

  wset, midplane_power_wid

  ;;-------------------------------------------------------
  ;; Now, begin menu selection
  ;;
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

     4: begin
        r_flag = 1-r_flag
        goto, plot_it
     end

     5: begin
        i_loglog = 1-i_loglog
        goto, plot_it
     end

     6: begin
        i_abs = 1-i_abs
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

     9: widget_control, midplane_power.top, /destroy

  endcase

  return

  plot_it:

  midplane_power_plot

  return
  
end
