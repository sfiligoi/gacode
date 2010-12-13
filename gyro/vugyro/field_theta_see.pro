pro field_theta_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_FIELD_THETA,$
    active_f,$
    widget

  active_f   = 1
  ;;-----------------------------------------------

  ;;----------------------------------------
  ;; first, check availablility
  ;;
  if exists_u eq 0 then return
  if xregistered('field_theta_see') then return
  ;;----------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  butt = widget_base(base,$
                     /row,$
                     /frame)

  w1 = widget_button(butt,$
                     value='All U',$
                     uvalue=1)

  w2 = widget_button(butt,$
                     value='U(t)',$
                     uvalue=2)

  w3 = widget_button(butt,         $
                     value='+t', $
                     uvalue=3)

  w4 = widget_button(butt,         $
                     value='-t', $
                     uvalue=4)

  w5 = widget_button(butt,         $
                     value='++t', $
                     uvalue=5)

  w6 = widget_button(butt,         $
                     value='--t', $
                     uvalue=6)

  w7 = widget_button(butt,         $
                     value='+n', $
                     uvalue=7)

  w8 = widget_button(butt,         $
                     value='-n', $
                     uvalue=8)

  w9 = widget_button(butt,         $
                     value='+r', $
                     uvalue=9)

  w10 = widget_button(butt,         $
                      value='-r', $
                      uvalue=10)

  w11 = widget_button(butt, $
                      value='field', $
                      uvalue=11)

  w12 = widget_button(butt, $
                      value='Done', $
                      uvalue=12)

  w13 = widget_button(butt, $
                      value='PS dump', $
                      uvalue=13)

  ;;-----------------------------------

  current_val = widget_base(base,$
                            /row,$
                            /frame)

  t_str = 't                 '
  t_label = widget_label(current_val,value=t_str)

  r_str = 'r                 '
  r_label = widget_label(current_val,value=r_str)

  n_str = 'n                 '
  n_label = widget_label(current_val,value=n_str)

  state = {t_label:t_label, $
           r_label:r_label, $
           n_label:n_label }

  ;;-----------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
      set_uvalue=state,$
      /no_copy, $
      /realize

  widget_control, draw, $
      get_value=widget

  xmanager,'field_theta_see', $
      base, $
      event='field_theta_event',$
      group_leader=group

end
