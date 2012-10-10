pro t_error_see, group=group

  common GLOBAL
  common T_ERROR_DATA

  ;;----------------------------------------
  if exists_t_error eq 0 then return    
  if xregistered('t_error_see') then return
  ;;----------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  row1 = widget_base(base,$
                     /row)

  x = widget_button(row1,$
                     value='Plot',$
                     uvalue=1)

  x = widget_button(row1,         $
                     value='+min', $
                     uvalue=2)

  x = widget_button(row1,         $
                     value='-min', $
                     uvalue=3)

  x = widget_button(row1,         $
                     value='PS dump', $
                     uvalue=4)

  x = widget_button(row1,         $
                     value='Done', $
                     uvalue=5)

  ;;-----------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

 ; widget_control, base,xoffset=500,yoffset=100, $
  widget_control, base, $
      set_uvalue=state,$
      /no_copy, $
      /realize

  widget_control, draw, $
      get_value=t_error_wid

  xmanager,'t_error_see', $
      base, $
      event='t_error_event',$
      group_leader=group

end
