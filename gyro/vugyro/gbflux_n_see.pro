pro gbflux_n_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_GBFLUX_N,$
    frac_flag,$
    bar_plot_flag,$
    widget

  frac_flag = 0
  bar_plot_flag = 0
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if (exists_gbflux_n eq 0) then return    
  if xregistered('gbflux_n_see') then return
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
                    uvalue=0)

  x = widget_button(row1, $
                    value='+field', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='+species', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='+moment', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='abs/frac', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='line/bar', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='units', $
                    uvalue=6)

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

  xmanager,'gbflux_n_see', $
    base, $
    event='gbflux_n_event',$
    group_leader=group

end
