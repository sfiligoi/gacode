pro dist_k_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_DIST_K,$
    widget
  ;;-----------------------------------------------

  ;;-------------------------------------
  ;; Return conditions
  ;;
  if exists_h eq 0 then return
  if xregistered('dist_k_see') then return
  ;;-------------------------------------

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
                    value='i/e', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='+ax', $
                    uvalue=7)

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
                    value='Done', $
                    uvalue=11)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
    set_uvalue=state, $
    /no_copy, $
    /realize

  widget_control, draw, $
    get_value=widget

  xmanager,'dist_k_see', $
    base, $
    event='dist_k_event',$
    group_leader=group

end
