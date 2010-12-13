pro midplane_n_power_event, midplane_n_power

  common GLOBAL
  common MIDPLANE_DATA

  widget_control, midplane_n_power.id, $
                  get_uvalue=uvalue

  wset, midplane_n_power_wid

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
        counter_up,in_c,n_n-1,1
        goto, plot_it
     end

     5: begin
        counter_dn,in_c,0,1 
        goto, plot_it
     end

     6: begin
        i_loglog = 1-i_loglog
        goto, plot_it
     end

     7: begin
        i_abs = 1-i_abs
        goto, plot_it
     end

     8: begin
        plot_mode = 2
        goto, plot_it
     end

     9: begin
        plot_export = 1
        goto, plot_it
     end

     10: widget_control, midplane_n_power.top, /destroy

  endcase

  return

  plot_it:

  midplane_n_power_plot

  return
  
end
