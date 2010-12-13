pro geometry_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_GEOMETRY,$
    i_geometry, $
    i_theta,$    
    widget

  i_geometry = 0
  i_theta = n_fine/2
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if (exists_geometry eq 0) then return
  if xregistered('geometry_see') then return
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
                    value='GEO FUNCTION',$
                    /menu)

  for i=0,n_geometry-1 do begin
     x1 = widget_button(x,$
                        value=geometry_label[i],$
                        uvalue=30+i)  
  endfor

  x = widget_button(row1, $
                    value='+i', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='-i', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=5)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base, $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
    set_uvalue=state,$
    /no_copy, $
    /realize

  widget_control, draw, $
    get_value=widget

  xmanager,'geometry_see', $
    base, $
    event='geometry_event',$
    group_leader=group

end
