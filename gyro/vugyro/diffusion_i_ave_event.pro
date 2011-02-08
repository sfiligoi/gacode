pro diffusion_i_ave_event, diffusion_i_ave

  common GLOBAL
  common PRIVATE_DIFFUSION_I_AVE  

  widget_control, diffusion_i_ave.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        i_spec = i_spec+1
        if (i_spec gt n_kinetic-1) then i_spec = 0
        goto, plot_it
     end   

     3: begin
        i_moment = 1-i_moment
        goto, plot_it
     end   

     4: begin
        i_f = i_f+1
        if (i_f ge n_field) then i_f = -1
        goto, plot_it
     end

     5: begin
        i_div = 1-i_div
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

     13: begin
        plot_export = 1
        goto, plot_it
     end

     14: begin
        plot_mode = 2
        goto, plot_it
     end

     15: widget_control, diffusion_i_ave.top, /destroy

  endcase
  
  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  diffusion_i_ave_plot

  return

end
