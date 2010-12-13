pro midplane_each_power_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_MIDPLANE_EACH_POWER,$
    bar_plot_n2_flag,$
    div_dn_flag,$
    widget

  bar_plot_n2_flag = 0
  div_dn_flag = 0
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if xregistered('midplane_each_power_see') then return
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
                    value='Plot', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='line/bar', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='div_del_n', $
                    uvalue=61)


  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='data', $
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
    get_value=widget

  xmanager,'midplane_each_power_see', $
    base, $
    event='midplane_each_power_event',$
    group_leader=group

end
