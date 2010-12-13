pro midplane_deriv_event, midplane_deriv

  common GLOBAL
  common MIDPLANE_DATA

  widget_control, midplane_deriv.id, $
                  get_uvalue=uvalue

  wset, midplane_deriv_wid

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
        counter_up,order,2,1
        goto, plot_it
     end

     5: begin
        counter_dn,order,0,1 
        goto, plot_it
     end

     6: begin
        log_flag = 1-log_flag
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

     9: widget_control, midplane_deriv.top, /destroy

  endcase

  return

  plot_it:

  midplane_deriv_plot

  return
  
end
