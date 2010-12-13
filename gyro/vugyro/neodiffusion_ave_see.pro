pro neodiffusion_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_NEODIFFUSION_AVE,$
    i_runave,$
    i_neomom,$
    zoom,$
    widget

  i_neomom = 0
  i_runave = 0
  zoom     = 1.0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_diff_neo eq 0 then return 
  if xregistered('diffusion_ave_see') then return
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
                    value='moment', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='run-ave', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='units', $
                    uvalue=6)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=7)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=8)

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

  xmanager,'neodiffusion_ave_see', $
           base, $
           event='neodiffusion_ave_event',$
           group_leader=group

end
