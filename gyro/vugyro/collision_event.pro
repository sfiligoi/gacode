pro collision_event, collision

  common GLOBAL
  common PLOT_VARIABLES
  common COLLISION_DATA

  widget_control, collision.id, $
                  get_uvalue=uvalue

  wset, collision_wid

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     1: goto, plot_it

     2: begin
        counter_up,i_sec,n_theta_section-1,1
        goto, plot_it
     end

     3: begin
        counter_dn,i_sec,0,1
        goto, plot_it
     end

     4: begin
        counter_up,ik_c,n_lambda-1,1
        goto, plot_it
     end

     5: begin 
        counter_dn,ik_c,0,1
        goto, plot_it
     end

     6: begin
        counter_up,i_c,n_r_write-1,1 
        goto, plot_it
     end

     7: begin
        counter_dn,i_c,0,1  
        goto, plot_it
     end

     8: begin
        plot_mode = 2
        goto, plot_it
     end

     9: widget_control, collision.top, /destroy
     
  endcase
  
  return

  plot_it:
  
  widget_control, collision.top, $
                  get_uvalue=state, $
                  /no_copy
  
  m_str = 'm = '+strtrim(string(i_sec),2)
  k_str = 'k = '+strtrim(string(ik_c),2)
  i_str = 'r = '+strtrim(string(i_c),2)

  widget_control, state.m_label, set_value=m_str
  widget_control, state.k_label, set_value=k_str
  widget_control, state.i_label, set_value=i_str      

  widget_control, collision.top, $
                  set_uvalue=state, $
                  /no_copy

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  plot_def_new,'collision'

  plot,[0],[0],/nodata,$
       xstyle=1,xrange=[-1,1],$
       ystyle=1,yrange=[-1,1],$
       xtitle='!4h!3',$
       ytitle='v!9!d#!3!n/v', $
       color=line

  oplot,xx[*],yy[*],psym=8,symsize=1.0

  for kp=0,n_lambda-1 do begin
     if mk[0,i_sec,kp,ik_c,0] gt 0 then begin

        oplot,[xx[mk[0,i_sec,kp,ik_c,0]-1,kp,0],$
               xx[mk[1,i_sec,kp,ik_c,0]-1,kp,0]],$
              [yy[mk[0,i_sec,kp,ik_c,0]-1,kp,0],$
               yy[mk[1,i_sec,kp,ik_c,0]-1,kp,0]],$
              psym=-8,symsize=2,color=color_vec[0]

     endif

     if mk[2,i_sec,kp,ik_c,0] gt 0 then begin
        oplot,[xx[mk[2,i_sec,kp,ik_c,0]-1,kp,0], $
               xx[mk[3,i_sec,kp,ik_c,0]-1,kp,0]],$
              [yy[mk[2,i_sec,kp,ik_c,0]-1,kp,0],$
               yy[mk[3,i_sec,kp,ik_c,0]-1,kp,0]],$
              psym=-8,symsize=2,color=color_vec[0]

     endif

  endfor

  oplot,[xx[i_sec,ik_c,0]],[yy[i_sec,ik_c,0]],$
        psym=8,$
        symsize=2,$
        color=color_vec[1]

  oplot,xx[i_sec,ik_c,0]*[1,1],[-1,1],$
        linestyle=line_vec[1]

  plot_finish

  return

end

