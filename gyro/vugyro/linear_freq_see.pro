pro linear_freq_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_LINEAR_FREQ, $
    widget, $
    real_flag, $
    zoom

  real_flag = 0
  zoom = 0.5
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions
  ;;
  if exists_omega eq 0 then return
  if xregistered('linear_freq_see') then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row)

  x = widget_button(row1,$
                    value='plot',$
                    uvalue=1)

  x = widget_button(row1,         $
                    value='+n', $
                    uvalue=2)

  x = widget_button(row1,         $
                    value='-n', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='Re(w)/Im(w)', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=7)

  x = widget_button(row1,         $
                    value='Done', $
                    uvalue=8)

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

  xmanager,'linear_freq_see', $
           base, $
           event='linear_freq_event',$
           group_leader=group

end
