pro neodiffusion_i_ave_event, neodiffusion_i_ave

  common GLOBAL
  common PRIVATE_NEODIFFUSION_I_AVE  

  widget_control, neodiffusion_i_ave.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     0: goto, plot_it

     1: begin
        i_spec = i_spec+1
        if (i_spec gt n_kinetic-1) then i_spec = 0
        goto, plot_it
     end   

     2: begin
        i_moment = 1+i_moment
        if (i_moment gt 2) then i_moment = 0
        goto, plot_it
     end   

     3: begin
        i_f = i_f+1
        if (i_f ge n_field) then i_f = -1
        goto, plot_it
     end

     20: begin
        i_div = 0
        goto, plot_it
     end

     21: begin
        i_div = 1
        goto, plot_it
     end

     22: begin
        i_div = 2
        goto, plot_it
     end

     30: begin
        i_turb_neo = 0
        goto, plot_it
     end

     31: begin
        i_turb_neo = 1
        goto, plot_it
     end

     32: begin
        i_turb_neo = 2
        goto, plot_it
     end

     33: begin
        i_turb_neo = 3
        goto, plot_it
     end

     34: begin
        i_turb_neo = 4
        goto, plot_it
     end

     35: begin
        i_turb_neo = 5
        goto, plot_it
     end

     5: begin
        i_zero = 1-i_zero
        goto, plot_it
     end
     
     9: begin
        zoom = zoom/2
        goto, plot_it
     end

     10: begin
        zoom = zoom*2
        goto, plot_it
     end

     11: begin
        i_units = 1-i_units
        goto, plot_it
     end
     
     12: begin
        n_ss = n_ss+dn_ss
        if (n_ss gt 4*dn_ss) then n_ss=0
        goto, plot_it 
     end

     13: begin
        n_ss = n_ss-dn_ss
        if (n_ss lt 0) then n_ss=4*dn_ss
        goto, plot_it 
     end

     14: begin
        plot_export = 1
        goto, plot_it
     end

     15: begin
        plot_mode = 2
        goto, plot_it
     end

     16: widget_control, neodiffusion_i_ave.top, /destroy

  endcase
  
  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  neodiffusion_i_ave_plot

  return

end
