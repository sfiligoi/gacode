pro diffusion_i_eff_event, diffusion_i_eff

  common GLOBAL
  common PRIVATE_DIFFUSION_I_EFF  

  widget_control, diffusion_i_eff.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     1: goto, plot_it

     4: begin
        i_f = i_f+1
        if (i_f ge n_field) then i_f = -1
        goto, plot_it
     end

     6: begin
        i_zero = 1-i_zero
        goto, plot_it
     end
     
     7: begin
        zoom = zoom/2
        goto, plot_it
     end

     8: begin
        zoom = zoom*2
        goto, plot_it
     end

     9: begin
        i_units = 1-i_units
        goto, plot_it
     end
     
     10: begin
        n_ss = n_ss+dn_ss
        if (n_ss gt 4*dn_ss) then n_ss=0
        goto, plot_it 
     end

     11: begin
        n_ss = n_ss-dn_ss
        if (n_ss lt 0) then n_ss=4*dn_ss
        goto, plot_it 
     end

     12: begin
        plot_export = 1
        goto, plot_it
     end

     13: begin
        plot_mode = 2
        goto, plot_it
     end

     14: widget_control, diffusion_i_eff.top, /destroy

  endcase
  
  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  diffusion_i_eff_plot

  return

end
