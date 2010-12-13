pro zerobar_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_ZEROBAR,$
    active_zerobar,$
    zoom,$
    widget

  active_zerobar = 1
  zoom = 1.0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_zerobar eq 0 then return    
  if xregistered('zerobar_see') then return
  ;;------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  butt = widget_base(base,$
                     /row,$
                     /frame)

  w1 = widget_button(butt, $
                     value='Function', $
                     /menu)

  wlist = widget_button(w1,value='Average phi',uvalue=1)  
  wlist = widget_button(w1,value='Average A',uvalue=2)  
  wlist = widget_button(w1,value='Average (d/dr) phi',uvalue=3)  
  wlist = widget_button(w1,value='Average rho (d/dr)^2 phi',uvalue=4)  

  w5 = widget_button(butt, $
                     value='ZOOM in', $
                     uvalue=5)

  w6 = widget_button(butt, $
                     value='ZOOM out', $
                     uvalue=6)

  w7 = widget_button(butt, $
                     value='PS dump', $
                     uvalue=7)

  w8 = widget_button(butt, $
                     value='Done', $
                     uvalue=8)

  ;;-----------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
      set_uvalue=state,$
      /no_copy, $
      /realize

  widget_control, draw, $
      get_value=widget

  xmanager,'zerobar_see', $
      base, $
      event='zerobar_event',$
      group_leader=group

end
