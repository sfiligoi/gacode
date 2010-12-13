pro zmoment_see, group=group

  common GLOBAL
  common ZMOMENT_DATA

  ;;------------------------------------------
  ;; Return conditions
  ;;    
  if exists_zmoment eq 0 then return    
  if xregistered('zmoment_see') then return
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
                    value='n/E', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='es/em', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=6)

  x = widget_button(row1,$
                    value='TYPE',$
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
                  get_value=zmoment_wid

  xmanager,'zmoment_see', $
           base, $
           event='zmoment_event',$
           group_leader=group

end
