pro osc_h_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_OSC_H,$
    gradient_flag,$
    pert_flag,$
    i_osc_mom,$
    zoom,$
    equil_flag,$
    widget

  pert_flag = 0
  i_osc_mom = 0
  gradient_flag = 0
  equil_flag = 1
  zoom = 1.0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_source eq 0 then return 
  if xregistered('osc_h_see') then return
  ;;------------------------------------------

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

  x = widget_button(row1,$
                    value='gradient',$
                    uvalue=2)

  x = widget_button(row1,$
                    value='moment',$
                    uvalue=3)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='EQ tog', $
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

  xmanager,'osc_h_see', $
    base, $
    event='osc_h_event',$
    group_leader=group

end
