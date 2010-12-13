pro linear_spec_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_LINEAR_SPEC,$
    widget, $
    real_flag,$
    zero_flag

  real_flag = 0
  zero_flag = 0
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions
  ;;
  if exists_omega eq 0 then return
  if xregistered('freq_spec_see') then return
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

  x = widget_button(row1, $
                    value='Re(w)/Im(w)', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='zero', $
                    uvalue=3)

  x = widget_button(row1,         $
                    value='PS dump', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=5)

  x = widget_button(row1,         $
                    value='Done', $
                    uvalue=6)

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

  xmanager,'linear_spec_see', $
    base, $
    event='linear_spec_event',$
    group_leader=group

end
