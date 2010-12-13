pro harmonics_see, group=group

  common GLOBAL

  ;;----------------------------------------
  ;; first, check availablility
  ;;
  if xregistered('harmonics_see') then return
  if exists_u eq 0 then return
  if simdir eq 'null' then return
  ;;----------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  butt = widget_base(base,$
                     /row,$
                     /frame)

  w1 = widget_button(butt,$
                     value='plot',$
                     uvalue=1)

  w2 = widget_button(butt,         $
                     value='+t', $
                     uvalue=2)

  w3 = widget_button(butt,         $
                     value='-t', $
                     uvalue=3)

  w4 = widget_button(butt,         $
                     value='++t', $
                     uvalue=4)

  w5 = widget_button(butt,         $
                     value='--t', $
                     uvalue=5)

  w6 = widget_button(butt,         $
                     value='+n', $
                     uvalue=6)

  w7 = widget_button(butt,         $
                     value='-n', $
                     uvalue=7)

  w8 = widget_button(butt,         $
                     value='+m', $
                     uvalue=8)

  w9 = widget_button(butt,         $
                     value='-m', $
                     uvalue=9)

  w10 = widget_button(butt, $
                      value='f tog', $
                      uvalue=10)

  w11 = widget_button(butt, $
                      value='PS dump', $
                      uvalue=11)

  w12 = widget_button(butt, $
                      value='export', $
                      uvalue=12)

  w13 = widget_button(butt, $
                      value='Done', $
                      uvalue=13)

  ;;-----------------------------------

  current_val = widget_base(base,$
                            /row,$
                            /frame)

  t_str = 't                 '
  t_label = widget_label(current_val,value=t_str)

  m_str = 'm                 '
  m_label = widget_label(current_val,value=m_str)

  n_str = 'n                 '
  n_label = widget_label(current_val,value=n_str)

  state = {t_label:t_label, $
           m_label:m_label, $
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
      get_value=harmonics_wid

  xmanager,'harmonics_see', $
      base, $
      event='harmonics_event',$
      group_leader=group

end
