pro field_r0_theta_see, group=group

  common GLOBAL

  ;;----------------------------------------
  ;; first, check availablility
  ;;
  if exists_field_r0 eq 0 then return
  if xregistered('field_r0_theta_see') then return
  ;;----------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row)

  x = widget_button(row1,$
                     value='Plot',$
                     uvalue=1)

  x = widget_button(row1,         $
                     value='+t', $
                     uvalue=2)

  x = widget_button(row1,         $
                     value='-t', $
                     uvalue=3)

  x = widget_button(row1,         $
                     value='++t', $
                     uvalue=4)

  x = widget_button(row1,         $
                     value='--t', $
                     uvalue=5)

  x = widget_button(row1,         $
                     value='+n', $
                     uvalue=6)

  x = widget_button(row1,         $
                     value='-n', $
                     uvalue=7)

  x = widget_button(row1, $
                      value='field', $
                      uvalue=8)

  x = widget_button(row1, $
                      value='PS dump', $
                      uvalue=9)

  x = widget_button(row1, $
                      value='Done', $
                      uvalue=10)

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
      get_value=field_r0_wid

  xmanager,'field_r0_theta_see', $
      base, $
      event='field_r0_theta_event',$
      group_leader=group

end
