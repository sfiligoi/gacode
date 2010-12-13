pro diffusion_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_DIFFUSION_AVE,$
    widget,$
    zoom, $
    i_tp

  zoom = 1.0
  i_tp = 0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if exists_diff eq 0 then return 
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
                    value='Plot', $
                    uvalue=0)

  x = widget_button(row1, $
                    value='i/e', $
                    uvalue=1)

  x = widget_button(row1, $
                    value='n/E', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='es/em', $
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
                    value='T+P', $
                    uvalue=7)

  x = widget_button(row1,$
                    value='TYPE',$
                    /menu)

  tlevels=['Line Plot','PDF']
  for i=0,1 do begin
     x1 = widget_button(x,$
                        value=tlevels[i],$
                        uvalue=20+i)  
  endfor

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=8)

  x = widget_button(row1, $
                    value='data', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=10)

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

  xmanager,'diffusion_ave_see', $
    base, $
    event='diffusion_ave_event',$
    group_leader=group

end
