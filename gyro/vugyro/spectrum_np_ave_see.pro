pro spectrum_np_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_SPECTRUM_NP_AVE,$
    widget
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if (exists_kxkyspec eq 0) then return    
  if xregistered("spectrum_np_ave_see") then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1,$
                    value='plot',$
                    uvalue=1)

  x = widget_button(row1,         $
                    value='+ax', $
                    uvalue=2)

  x = widget_button(row1,         $
                    value='-ax', $
                    uvalue=3)

  x = widget_button(row1,         $
                    value='+az', $
                    uvalue=4)

  x = widget_button(row1,         $
                    value='-az', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='p^2', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='log-tog', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=8)

  x = widget_button(row1,         $
                    value='Done', $
                    uvalue=9)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
                  /no_copy, $
                  /realize

  widget_control, draw, $
                  get_value=widget

  xmanager,'spectrum_np_ave_see', $
           base, $
           event='spectrum_np_ave_event',$
           group_leader=group

end
