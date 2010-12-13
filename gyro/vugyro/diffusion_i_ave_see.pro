pro diffusion_i_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_DIFFUSION_I_AVE,$
    i_div,$
    neo_flag,$
    zoom,$
    widget

  i_div    = 0
  zoom     = 5.0
  neo_flag = 0
  ;;-----------------------------------------------

  ;;-----------------------------------------------
  ;; Return conditions
  ;;
  if exists_diff_i eq 0 then return 
  if xregistered('diffusion_i_ave_see') then return
  ;;-----------------------------------------------

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
                    value='i/e', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='n/E', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='es/em', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='divisor', $
                    uvalue=5)
  
  x = widget_button(row1, $
                    value='zero', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=8)


  x =  widget_button(row1, $
                     value='units', $
                     uvalue=9)

  row2 = widget_base(base, $
                     /row, $
                     /frame)

  x = widget_button(row2, $
                    value='ss+', $
                    uvalue=10)

  x = widget_button(row2, $
                    value='ss-', $
                    uvalue=11)

  x = widget_button(row2, $
                    value='add neo', $
                    uvalue=12 )

  x = widget_button(row2, $
                    value='data', $
                    uvalue=13)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=14)

  x = widget_button(row2, $
                    value='Done', $
                    uvalue=15)

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

  xmanager,'diffusion_i_ave_see', $
    base, $
    event='diffusion_i_ave_event',$
    group_leader=group

end
