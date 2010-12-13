pro dist_k_event, dist_k

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_DIST_K
  
  widget_control, dist_k.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it
     
     2: begin
        counter_up,t_c,n_time1,1
        goto, plot_it
     end

     3: begin
        counter_dn,t_c,0,1
        goto, plot_it
     end

     4: begin
        counter_up,t_c,n_time1,10 
        goto, plot_it
     end

     5: begin 
        counter_dn,t_c,0,10
        goto, plot_it
     end

     6: begin
        i_spec = i_spec+1
        if i_spec gt n_kinetic-1 then i_spec = 0
        goto, plot_it
     end
     
     7: begin
        counter_up,ax_c,90,5
        goto, plot_it
     end

     8: begin
        counter_dn,ax_c,-90,5
        goto, plot_it
     end

     9: begin
        counter_up,az_c,90,5
        goto, plot_it
     end

     10: begin 
        counter_dn,az_c,-90,5
        goto, plot_it
     end

     11: widget_control, dist_k.top, /destroy

  endcase

  return

  plot_it:
  
  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  zmax = 1.1*max(abs(ht[*,*,*,*,t_c]))

  ymin = 0.0
  ymax = max(lambda)

  label = 'null'
  plot_def_new,label

  surface,hp[0,*,*,i_spec,t_c],theta_t_p/!pi,xi_p,$
    zstyle=1,zrange=[-1,1]*zmax,$
    xstyle=1,xrange=[-1,1],$
    ystyle=1,yrange=[-1,1],$
    ax=ax_c,az=az_c,title="Hp t="+strtrim(string(t[t_c]),2), $
    ytitle='k',xtitle='!4h/2p!3',color=color_vec[0]

  surface,hp[1,*,*,i_spec,t_c],theta_t_p/!pi,-xi_p,$
    zstyle=1,zrange=[-1,1]*zmax,$
    xstyle=1,xrange=[-1,1],$
    ystyle=1,yrange=[-1,1],$
    ax=ax_c,az=az_c,title="Hp t="+strtrim(string(t[t_c]),2), $
    ytitle='k',xtitle='!4h/2p!3',color=color_vec[1],$
    /noerase

  surface,ht[0,*,*,i_spec,t_c],theta_t_t/!pi,xi_t,$
    zstyle=1,zrange=[-1,1]*zmax,$
    xstyle=1,xrange=[-1,1],$
    ystyle=1,yrange=[-1,1],$
    ax=ax_c,az=az_c,title="Hp t="+strtrim(string(t[t_c]),2), $
    ytitle='k',xtitle='!4h/2p!3',color=color_vec[0],$
    /noerase

  surface,ht[1,*,*,i_spec,t_c],theta_t_t/!pi,-xi_t,$
    zstyle=1,zrange=[-1,1]*zmax,$
    xstyle=1,xrange=[-1,1],$
    ystyle=1,yrange=[-1,1],$
    ax=ax_c,az=az_c,title="Hp t="+strtrim(string(t[t_c]),2), $
    ytitle='k',xtitle='!4h/2p!3',color=color_vec[1],$
    /noerase

  plot_finish

  return

end
