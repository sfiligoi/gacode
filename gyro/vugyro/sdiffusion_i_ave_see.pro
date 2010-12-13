pro sdiffusion_i_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_SDIFFUSION_I_AVE,$
    zoom,$
    widget, $
    i_tp

  zoom = 5.0
  i_tp = 0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_sdiff_i eq 0 then return 
  if xregistered('sdiffusion_i_ave_see') then return
  ;;------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='i/e', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='M/H', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='units', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='T+P', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='zero', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=9)

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

  xmanager,'sdiffusion_i_ave_see', $
           base, $
           event='sdiffusion_i_ave_event',$
           group_leader=group

end
