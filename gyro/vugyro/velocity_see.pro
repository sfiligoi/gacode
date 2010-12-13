pro velocity_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_VELOCITY,$
    widget,$
    i_plot

  i_plot = 0
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if exists_velocity eq 0 then return    
  if xregistered('velocity_see') then return
  ;;----------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row)

  x = widget_button(row1,$
                    value='Plot',$
                    uvalue=1)

  x = widget_button(row1,         $
                    value='+t', $
                    uvalue=2)

  x = widget_button(row1,         $
                    value='-t', $
                    uvalue=3)

  x = widget_button(row1,         $
                    value='++t', $
                    uvalue=4)

  x = widget_button(row1,         $
                    value='--t', $
                    uvalue=5)

  x = widget_button(row1,         $
                    value='+n', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='-n', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='i/e', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='n/E', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='es/em', $
                    uvalue=10)

  x = widget_button(row1,$
                    value='TYPE',$
                    /menu)

  tlevels=['D(e,l)','D(e)','D(l)']
  for i=0,2 do begin
     x1 = widget_button(x,$
                        value=tlevels[i],$
                        uvalue=20+i)  
  endfor

  x = widget_button(row1,         $
                    value='PS dump', $
                    uvalue=11)

  x = widget_button(row1,         $
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

  xmanager,'velocity_see', $
    base, $
    event='velocity_event',$
    group_leader=group

end
