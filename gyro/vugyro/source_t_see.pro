pro source_t_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_SOURCE_T,$
    zoom,$
    widget

  zoom = 1.0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_source eq 0 then return    
  if xregistered('source_t_see') then return
  ;;------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1,$
                    value='Plot',$
                    uvalue=1)

  x = widget_button(row1,$
                    value='spec',$
                    uvalue=2)  

  x = widget_button(row1,$
                    value='moment',$
                    uvalue=3)  

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=7)

  x = widget_button(row1, $
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

  xmanager,'source_t_see', $
    base, $
    event='source_t_event',$
    group_leader=group

end
