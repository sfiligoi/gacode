pro each_diffusion_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_EACH_DIFFUSION,$
    frac_flag,$
    bar_plot_flag,$
    i_QL_flag,$
    i_QL_mult,$
    i_QL_g_phi,$
    widget

  frac_flag = 0
  bar_plot_flag = 0
  i_QL_flag = 0
  i_QL_mult = 0
  i_QL_g_phi = 0
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if (exists_diff_n eq 0) then return    
  if xregistered('each_diffusion_see') then return
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
                    value='README', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='QL', $
                    uvalue=11)

  x = widget_button(row1, $
                    value='QL-mult', $
                    uvalue=12)

  x = widget_button(row1, $
                    value='QL-norm', $
                    uvalue=13)

  x = widget_button(row1, $
                    value='i/e', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='n/E', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='es/em', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='abs/frac', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='line/bar', $
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

  xmanager,'each_diffusion_see', $
    base, $
    event='each_diffusion_event',$
    group_leader=group

end
