pro osc_see, group=group

  common GLOBAL
  common ZMOMENT_DATA

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_OSC,$
    i_gradient,$
    i_ss_mn_den_plot,$
    i_ss_mn,$
    n_ss_bin,$
    zoom,$
    i_div,$
    equil_flag,$
    widget

  i_div      = 0
  i_ss_mn_den_plot = 0
  i_ss_mn    = 0
  i_gradient = 0
  n_ss_bin   = 20
  equil_flag     = 1
  zoom       = 1.0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_zmoment eq 0 then return    
  if xregistered('osc_see') then return
  ;;------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  ;; ROW 1

  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='i/e', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='moment', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='gradient', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='zero', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='EQ tog', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='ss+', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='ss-', $
                    uvalue=11)

  ;; ROW 2

  row2 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row2, $
                    value='ss density tog', $
                    uvalue=12)

  x = widget_button(row2, $
                    value='bin++', $
                    uvalue=13)

  x = widget_button(row2, $
                    value='bin--', $
                    uvalue=14)

  x = widget_button(row2, $
                    value='surf/mode tog', $
                    uvalue=15)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=16)

  x = widget_button(row2, $
                    value='data', $
                    uvalue=17)

  x = widget_button(row2, $
                    value='Done', $
                    uvalue=18)

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

  xmanager,'osc_see', $
    base, $
    event='osc_event',$
    group_leader=group

end
