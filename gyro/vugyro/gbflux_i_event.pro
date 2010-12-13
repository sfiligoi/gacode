pro gbflux_i_event, gbflux_i_wid

  common GLOBAL
  common PRIVATE_GBFLUX_I  

  widget_control, gbflux_i_wid.id, get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        i_spec = i_spec+1
        if (i_spec ge n_kinetic) then i_spec = 0
        goto, plot_it
     end   

     3: begin
        i_moment = i_moment+1
        if (i_moment ge p_moment) then i_moment = 0
        goto, plot_it
     end   

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
        if (exists_exp_derived eq 1) then begin
           i_units = i_units+1
           if (i_units eq 3) then i_units = 0
        endif
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

     15: widget_control, gbflux_i_wid.top, /destroy

  endcase
  
  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  gbflux_i_plot

  return

end
