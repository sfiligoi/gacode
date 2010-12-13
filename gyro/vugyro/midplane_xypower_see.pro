pro midplane_xypower_see, group=group

  common GLOBAL
  common MIDPLANE_DATA
  common POLOIDAL_DATA

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
  print, 'n_pwr=',n_pwr
  if n_n eq 1 then return
  if n_pwr eq 0 then return    
  if xregistered('power_spectrum_see') then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1, $
                    value='plot', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='index up', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='index dn', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='quad/lin', $
                    uvalue=4)


  x = widget_button(row1, $
                    value='lin/log', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=6)

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
                  get_value=midplane_xypower_wid

  xmanager,'midplane_xypower_see', $
           base, $
           event='midplane_xypower_event',$
           group_leader=group

end
