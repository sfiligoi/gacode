pro diffusion_i_eff_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_DIFFUSION_I_EFF,$
    i_div,$
    zoom,$
    widget

  i_div  = 0
  zoom   = 5.0
  ;;-----------------------------------------------

  ;;-----------------------------------------------
  ;; Return conditions
  ;;
  if exists_diff_i eq 0 then return 
  if n_kinetic eq n_ion then return
  if xregistered('diffusion_i_eff_see') then return
  ;;-----------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row)

  x = widget_button(row1, $
                    value='Plot', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='es/em', $
                    uvalue=4)

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

  x = widget_button(row1, $
                    value='ss+', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='ss-', $
                    uvalue=11)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=12)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=13)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=14)

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

  xmanager,'diffusion_i_eff_see', $
    base, $
    event='diffusion_i_eff_event',$
    group_leader=group

end
