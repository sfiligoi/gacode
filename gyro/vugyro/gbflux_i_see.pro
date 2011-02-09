pro gbflux_i_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_GBFLUX_I,$
    i_div,$
    zoom,$
    widget

  i_div    = 0
  zoom     = 1.0
  ;;-----------------------------------------------

  ;;-----------------------------------------------
  ;; Return conditions
  ;;
  if exists_gbflux_i eq 0 then return 
  if xregistered('diffusion_i_ave_see') then return
  ;;-----------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='+field', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='+species', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='+moment', $
                    uvalue=3)
  
  x = widget_button(row1, $
                    value='zero', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=8)

  x =  widget_button(row1, $
                     value='units', $
                     uvalue=9)

  row2 = widget_base(base, $
                     /row, $
                     /frame)

  x = widget_button(row2, $
                    value='ss+', $
                    uvalue=10)

  x = widget_button(row2, $
                    value='ss-', $
                    uvalue=11)

  x = widget_button(row2, $
                    value='data', $
                    uvalue=13)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=14)

  x = widget_button(row2, $
                    value='Done', $
                    uvalue=15)

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

  xmanager,'gbflux_i_see', $
    base, $
    event='gbflux_i_event',$
    group_leader=group

end
