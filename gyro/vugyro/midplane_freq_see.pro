pro midplane_freq_see, group=group

  common GLOBAL
  common MIDPLANE_DATA
  common POLOIDAL_DATA

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common  PRIVATE_MIDPLANE_FREQ,$
    spec_amp,$
    widget

  spec_amp = 1.0
  ;;-----------------------------------------------

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
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
                    value='lin/log', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='+zoom', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='-zoom', $
                    uvalue=7)

  x = widget_button(row1,$
                    value='TYPE',$
                    uvalue=8)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=11)

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

  xmanager,'midplane_freq_see', $
    base, $
    event='midplane_freq_event',$
    group_leader=group

end
