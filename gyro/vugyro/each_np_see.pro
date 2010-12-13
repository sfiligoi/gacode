pro each_np_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_EACH_NP,$
    zoom,$
    widget,$
    t_dot, $
    p_c

  zoom  = 1.0
  t_dot = 0
  p_c   = 1
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if (exists_kxkyspec eq 0) then return    
  if xregistered('each_np_see') then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1,$
                    value='plot',$
                    uvalue=1)

  x = widget_button(row1,$
                    value='+n',$
                    uvalue=2)

  x = widget_button(row1,$
                    value='-n',$
                    uvalue=3)

  x = widget_button(row1,$
                    value='+p',$
                    uvalue=4)

  x = widget_button(row1,$
                    value='-p',$
                    uvalue=5)

  x = widget_button(row1,$
                    value='p^2',$
                    uvalue=6)

  x = cw_field(row1,$
               /return_events, $
               title='ZOOM', $
               /floating, $
               value=1.0, $
               xsize=8, $
               uvalue=7)

  x = widget_button(row1,         $
                    value='log-tog', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='+dot', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='-dot', $
                    uvalue=10)

  row2 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=11)

  x = widget_button(row2,         $
                    value='Done', $
                    uvalue=12)

  n_str = ' n      '
  n_label = widget_label(row2,value=n_str)
  p_str = ' p      '
  p_label = widget_label(row2,value=p_str)

  state = {n_label:n_label, $
           p_label:p_label}

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
                  set_uvalue=state,$
                  /no_copy, $
                  /realize

  widget_control, draw, $
                  get_value=widget

  xmanager,'each_np_see', $
           base, $
           event='each_np_event',$
           group_leader=group

end
