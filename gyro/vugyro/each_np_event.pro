pro each_np_event, each_np

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_EACH_NP
  
  widget_control, each_np.id, $
    get_uvalue=uvalue

  wset, widget
  
  ;;-------------------------------------------------------
  ;; Now, begin menu selection
  ;;
  case (uvalue) of 
     
     1: goto, plot_it   

     2: begin
        counter_up,in_c,n_n-1,1
        goto, plot_it
     end

     3: begin
        counter_dn,in_c,0,1
        goto, plot_it
     end

     4: begin
        counter_up,p_c,n_r-1,1
        goto, plot_it
     end

     5: begin
        counter_dn,p_c,0,1
        goto, plot_it
     end

     6: begin
        p_squ=1-p_squ
        goto, plot_it
     end

     7: begin
        zoom = each_np.value
        goto, plot_it
     end

     8: begin 
        i_log = 1-i_log
        goto, plot_it
     end

     9: begin
        counter_up,t_dot,n_time1,1 
        goto, plot_it
     end

     10: begin
        counter_dn,t_dot,0,1
        goto, plot_it
     end

     11: begin
        plot_mode = 2
        goto, plot_it
     end

     12: widget_control, each_np.top, /destroy

  endcase                              

  ;;----------------------------------------------------------

  return

  plot_it:

  str_p = strtrim(string(fix(p_c-n_r/2)),2)
  str_n = strtrim(string(fix(n_tor[in_c])),2)

  plot_def_new,'each_np_'+str_p+'_'+str_n

  title  = '!4u!3(n='+str_n+',p='+str_p+')!3'

  ;;-------------------------------------------------
  ;; widget control for label functions
  ;;
  widget_control, each_np.top, $
    get_uvalue=state, $
    /no_copy
  
  n_str = ' n = '+str_n
  p_str = ' p = '+str_p

  widget_control, state.n_label, set_value=n_str  
  widget_control, state.p_label, set_value=p_str  

  widget_control, each_np.top, $
    set_uvalue=state, $
    /no_copy
  ;;-------------------------------------------------

  y = fltarr(n_time)
  y[*] = sqrt(kxkyspec[p_c,in_c,*])

  if (p_squ eq 1) then begin
     p_number = p_c-n_r/2
     print, 'p_number=',p_number
     print, 'p_squ=',p_squ

     y[*] = y[*]*p_number*p_number
     if (p_number eq 0) then begin
        y[*]=0.0
        for p=0,n_r-1 do begin
           y[*]=y[*]+kxkyspec[p,in_c,*]*(p-n_r/2)^4
        endfor
        y[*] = sqrt(y[*])/n_r
     endif
  endif

  ymax = max(y[1:n_time1])
  ymin = min(y[1:n_time1])

  if i_log eq 1 then begin

     plot_io,[0],[1],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(t),max(t)],$
       xtitle=csa_string,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax]/zoom,$
       color=line

  endif else begin

     plot,[1],[1],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(t),max(t)],$
       xtitle=csa_string,$
       ystyle=1,$
       yminor=0,$
       yrange=[0,ymax]/zoom,$
       color=line

  endelse

  oplot,t,y,color=color_vec[0]   

  oplot,[1,1]*t[t_dot],[1,1]*y[t_dot],$
    color=line,$
    psym=-8,$
    symsize=dotsize

  oplot,[t[0],t[n_time1]],[1,1]*y[t_dot],$
    color=line,$
    linestyle=1,$
    symsize=dotsize

  plot_finish

  return

end
