pro flow_i_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_FLOW_I_AVE,$
    active_flow_i_ave,$
    neo_flag,$
    widget,$
    zoom

  active_flow_i_ave = 3
  zoom              = 10.0
  neo_flag          = 0
  ;;-----------------------------------------------

  ;;-----------------------------------------------
  ;; Return conditions
  ;;
  if exists_diff_i eq 0 then return    
  if xregistered('flow_i_ave_see') then return
  ;;-----------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  x = widget_button(row1,$
                    value='Function',$
                    /menu)

  x1 = widget_button(x,value='Ion energy',uvalue=1)  
  x1 = widget_button(x,value='Electron energy',uvalue=2)  
  x1 = widget_button(x,value='Plasma flow',uvalue=3)
  x1 = widget_button(x,value='Momentum flow',uvalue=4)
  x1 = widget_button(x,value='Ion heating flow',uvalue=5)
  x1 = widget_button(x,value='Electron heating flow',uvalue=6)
  x1 = widget_button(x,value='Total heating flow',uvalue=7)

  x = widget_button(row1, $
                    value='add neo', $
                    uvalue=111)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='zero', $
                    uvalue=91)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=11)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=12)

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

  xmanager,'flow_i_ave_see', $
    base, $
    event='flow_i_ave_event',$
    group_leader=group

end
