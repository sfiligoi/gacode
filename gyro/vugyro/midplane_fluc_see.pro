pro midplane_fluc_see, group=group

  common GLOBAL
  common MIDPLANE_DATA
  common POLOIDAL_DATA

  ;;----------------------------------------------------------
  ;; Return conditions
  ;;
  if n_n eq 1 then return
  if n_pwr eq 0 then return
  if xregistered('midplane_fluc_see') then return
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

  x = widget_button(row1,         $
                    value='+t', $
                    uvalue=4)

  x = widget_button(row1,         $
                    value='-t', $
                    uvalue=5)

  x = widget_button(row1,         $
                    value='++t', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='--t', $
                    uvalue=7)

  x = widget_button(row1, $
                      value='refine', $
                      /menu)

  for i=0,n_mUt-1 do begin
    x1 = widget_button(x,$
                          value=s_mUt[i],$
                          uvalue=40+i)  
  endfor

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='movie', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=10)

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
                  get_value=midplane_fluc_wid

  xmanager,'midplane_fluc_see', $
           base, $
           event='midplane_fluc_event',$
           group_leader=group

end
