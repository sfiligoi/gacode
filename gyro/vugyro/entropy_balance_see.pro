pro entropy_balance_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_ENTROPY_BALANCE,$
    entropy_balance_max,$
    widget

  entropy_balance_max = max(entropy)
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if exists_entropy eq 0 then return    
  if xregistered('entropy_balance_see') then return
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

  x = widget_button(row1, $
                    value='+species', $
                    uvalue=2)

  x = widget_button(row1,         $
                    value='+min', $
                    uvalue=3)

  x = widget_button(row1,         $
                    value='-min', $
                    uvalue=4)

  x = widget_button(row1,         $
                    value='PS dump', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='Done', $
                    uvalue=7)

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

  xmanager,'entropy_balance_see', $
           base, $
           event='entropy_balance_event',$
           group_leader=group

end
