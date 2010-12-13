pro spectrum_np_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_SPECTRUM_NP,$
    widget
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if (exists_kxkyspec eq 0) then return    
  if xregistered("spectrum_np_see") then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  ;; ROW 1

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
                    value='+ax', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='-ax', $
                    uvalue=8)

  x = widget_button(row1,         $
                    value='+az', $
                    uvalue=9)

  x = widget_button(row1,         $
                    value='-az', $
                    uvalue=10)

  x = widget_button(row1,         $
                    value='p^2', $
                    uvalue=7)

  x = widget_button(row1,         $
                    value='log-tog', $
                    uvalue=11)

  ;; ROW 2

  row2 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=12)

  x = widget_button(row2, $
                    value='Done', $
                    uvalue=13)

  t_str = ' (null)           '
  t_label = widget_label(row2,value=t_str)

  state = {t_label:t_label}

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

  xmanager,'spectrum_np_see', $
           base, $
           event='spectrum_np_event',$
           group_leader=group

end
