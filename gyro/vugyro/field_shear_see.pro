;; setup: plot radial shear in (phi,A)

pro field_shear_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_FIELD_SHEAR,$
    active_f,$
    widget

  active_f   = 1
  ;;-----------------------------------------------

  ;;-------------------------------------------
  ;; Return conditions
  ;;
  if exists_u eq 0 then return
  if xregistered('field_shear_see') then return
  ;;-------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  butt = widget_base(base,$
                     /row)

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
                     value='+j', $
                     uvalue=9)

  w10 = widget_button(butt,         $
                      value='-j', $
                      uvalue=10)

  ;;-----------------------------------

  current_val = widget_base(base,$
                            /row, $
                            /frame)

  w13 = widget_button(butt, $
                      value='field', $
                      uvalue=13)

  w14 = widget_button(current_val,         $
                      value='Done', $
                      uvalue=14)
  
  w15 = widget_button(current_val,         $
                      value='PS dump', $
                      uvalue=15)

  t_str = 't                  '
  t_label = widget_label(current_val,value=t_str)

  j_str = 'theta/pi           '
  j_label = widget_label(current_val,value=j_str)

  n_str = 'n                  '
  n_label = widget_label(current_val,value=n_str)

  state = {t_label:t_label, $
           j_label:j_label, $
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

  xmanager,'field_shear_see', $
    base, $
    event='field_shear_event',$
    group_leader=group

end
