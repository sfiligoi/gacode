pro osc_event, osc

  common GLOBAL
  common ZMOMENT_DATA
  common PRIVATE_OSC
  
  widget_control, osc.id, $
    get_uvalue=uvalue

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
;  need to restrict the range of i_moment to 0-1:
        i_moment = 1 < ( 1-i_moment > 0 )
        goto, plot_it
     end

     4: begin
        i_gradient = 1-i_gradient
        goto, plot_it
     end
     
     6: begin
        zoom = zoom/2
        goto, plot_it
     end

     7: begin
        zoom = zoom*2
        goto, plot_it
     end

     8: begin
        i_zero = 1-i_zero
        goto, plot_it
     end


     9: begin
        equil_flag = 1-equil_flag
        goto, plot_it
     end

     10: begin
        n_ss = n_ss+dn_ss
        if (n_ss gt (n_n-1)*dn_ss) then n_ss=0
        goto, plot_it 
     end

     11: begin
        n_ss = n_ss-dn_ss
        if (n_ss lt 0) then n_ss=(n_n-1)*dn_ss
        goto, plot_it 
     end

     12: begin
        i_ss_mn_den_plot = 1-i_ss_mn_den_plot
        goto, plot_it
     end

     13: begin
        n_ss_bin = n_ss_bin*2
        goto, plot_it
     end

     14: begin
        n_ss_bin = n_ss_bin/2
        if (n_ss_bin lt 0) then n_ss_bin=1
        goto, plot_it
     end

     15: begin
        i_ss_mn = 1-i_ss_mn 
        goto, plot_it
     end

     16: begin
        plot_mode = 2
        goto, plot_it
     end

     17: begin
        plot_export = 1
        goto, plot_it
     end

     18: widget_control, osc.top, /destroy

  endcase
  ;;------------------------------------------------------------------

  return

  plot_it:

  osc_plot

  return

end
