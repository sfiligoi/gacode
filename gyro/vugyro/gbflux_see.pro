pro gbflux_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_GBFLUX,$
    widget,$
    zoom, $
    i_tp

  zoom = 1.0
  i_tp = 0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_gbflux eq 0 then return 
  if xregistered('gbflux_see') then return
  ;;------------------------------------------

  base = widget_base(title=simdir,/column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=0)

  x = widget_button(row1, $
                    value='+field', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='+species', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='+moment', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='trap/pass', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='units', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=10)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,xsize=sx,ysize=sy)

  widget_control, base, $
    set_uvalue=state,$
    /no_copy, $
    /realize

  widget_control, draw, $
    get_value=widget

  xmanager,'gbflux_see', $
    base, $
    event='gbflux_event',$
    group_leader=group

end
