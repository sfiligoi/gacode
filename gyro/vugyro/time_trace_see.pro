pro time_trace_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_TIME_TRACE,$
    i_er_shear,$
    active_f,$
    t_dot,$
    index,$
    widget

  index = 1
  i_er_shear = 0
  active_f   = 1
  t_dot = 0
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if exists_u eq 0 then return
  if xregistered('time_trace_see') then return
  ;;----------------------------------------------------------

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
                    value='U(i,j)',$
                    uvalue=1)

  x = widget_button(row1,$
                    value='U(i)',$
                    uvalue=2)

  x = widget_button(row1,$
                    value='U(j)',$
                    uvalue=3)

  x = widget_button(row1,$
                    value='|U|',$
                    uvalue=4)

  x = widget_button(row1,         $
                    value='+j', $
                    uvalue=5)

  x = widget_button(row1,         $
                    value='-j', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='+i', $
                    uvalue=7)

  x = widget_button(row1,         $
                    value='-i', $
                    uvalue=8)

  x = widget_button(row1,$
                    value='Er_shear',$
                    uvalue=81)

  x = widget_button(row1,         $
                    value='+n', $
                    uvalue=9)

  x = widget_button(row1,         $
                    value='-n', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='log/lin', $
                    uvalue=11)

  x = widget_button(row1, $
                    value='+dot', $
                    uvalue=12)

  x = widget_button(row1, $
                    value='-dot', $
                    uvalue=13)

  ;; Row 2

  row2 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row2, $
                    value='field', $
                    uvalue=14)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=15)

  x = widget_button(row2, $
                    value='data', $
                    uvalue=16)

  x = widget_button(row2, $
                    value='Done', $
                    uvalue=17)

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

  xmanager,'time_trace_see', $
    base, $
    event='time_trace_event',$
    group_leader=group

end
