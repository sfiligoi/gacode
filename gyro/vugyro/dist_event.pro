pro dist_event, dist

  common GLOBAL
  common PRIVATE_DIST
  
  widget_control, dist.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it
     
     2: begin
        t_c = t_c+1
        if t_c gt n_time1 then t_c = n_time1 
        goto, plot_it
     end

     3: begin
        t_c = t_c-1
        if t_c lt 0 then t_c = 0 
        goto, plot_it
     end

     4: begin
        t_c = t_c+10
        if t_c gt n_time1 then t_c = n_time1 
        goto, plot_it
     end

     5: begin
        t_c = t_c-10
        if t_c lt 0 then t_c = 0 
        goto, plot_it
     end

     6: begin
        ik_c = ik_c+1
        if ik_c gt n_lambda-1 then ik_c = n_lambda-1 
        goto, plot_it
     end

     7: begin
        ik_c = ik_c-1
        if ik_c lt 0 then ik_c = 0
        goto, plot_it
     end

     8: begin
        i_spec = i_spec+1
        if (i_spec ge n_kinetic) then i_spec = 0
        goto, plot_it
     end

     9: begin
        dist_axis = 1-dist_axis
        goto, plot_it
     end

     10: begin
        plot_mode = 2
        goto, plot_it
     end

     11: widget_control, dist.top, /destroy
     
  endcase

  return

  plot_it:
  
  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  dist_plot

  return

end
