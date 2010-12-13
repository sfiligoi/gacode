pro balloon_see, group=group

  common GLOBAL
  common PRIVATE_BALLOON

  ;;------------------------------------------;
  ;; Return conditions;
  ;;
  if (exists_balloon eq 0) then return
  if xregistered('balloon_see') then return
  ;;----------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  ;; Row 1

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1,$
                    value='plot',$
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

  x = widget_button(row1, $
                    value='+field', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='+l0', $
                    uvalue=7)

  x = widget_button(row1,           $
                    value='+norm',  $
                    uvalue=8)

  x = widget_button(row1,          $
                    value='-norm',  $
                    uvalue=9)

  x = widget_button(row1,           $
                    value='+cells',  $
                    uvalue=10)

  x = widget_button(row1,          $
                    value='-cells',  $
                    uvalue=11)

  x = widget_button(row1, $
                    value='+vert', $
                    uvalue=12)

  x = widget_button(row1, $
                    value='-vert', $
                    uvalue=13)

  ;; Row 2

  row2 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=14)

  x = widget_button(row2, $
                    value='export', $
                    uvalue=15)

  x = widget_button(row2, $
                    value='MPEG', $
                    uvalue=16)

  x = widget_button(row2, $
                    value='Done', $
                    uvalue=17)

  t_str = '  (null)      '
  t_label = widget_label(row2,value=t_str)

  n_str = '  (null)                      '
  n_label = widget_label(row2,value=n_str)

  state = { t_label:t_label , n_label:n_label }

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

  xmanager,'balloon_see', $
    base, $
    event='balloon_event',$
    group_leader=group

end
