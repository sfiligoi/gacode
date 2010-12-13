pro field_r_see, group=group

  common GLOBAL

  ;;----------------------------------------------------------
  ;; Return conditions
  ;;
  if exists_u eq 0 then return
  if xregistered('field_r_see') then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  ;; Row 1

  row1 = widget_base(base,$
                     /row)
  
  x = widget_button(row1,$
                    value='Plot',$
                    uvalue=1)

  x = widget_button(row1,$
                    value='one/all',$
                    uvalue=2)

  x = widget_button(row1,         $
                    value='+t', $
                    uvalue=3)

  x = widget_button(row1,         $
                    value='-t', $
                    uvalue=4)

  x = widget_button(row1,         $
                    value='++t', $
                    uvalue=5)

  x = widget_button(row1,         $
                    value='--t', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='+n', $
                    uvalue=7)

  x = widget_button(row1,         $
                    value='-n', $
                    uvalue=8)

  x = widget_button(row1,         $
                    value='+j', $
                    uvalue=9)

  x = widget_button(row1,         $
                    value='-j', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='field', $
                    uvalue=11)

  x = widget_button(row1, $
                    value='|U|', $
                    uvalue=12)

  x = widget_button(row1, $
                    value='log10_|U|', $
                    uvalue=13)

  ;; Row 2

  row2 = widget_base(base,$
                     /row, $
                     /frame)

  x = widget_button(row2,         $
                    value='PS dump', $
                    uvalue=14)

  x = widget_button(row2,         $
                    value='export', $
                    uvalue=15)

  x = widget_button(row2,         $
                    value='Done', $
                    uvalue=16)
  

  j_str = ' (null)               '
  j_label = widget_label(row2,value=j_str)

  state = {j_label:j_label}

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
                  get_value=field_r_wid

  xmanager,'field_r_see', $
           base, $
           event='field_r_event',$
           group_leader=group

end
