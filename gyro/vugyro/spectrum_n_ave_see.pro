pro spectrum_n_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_SPECTRUM_N_ave,$
    zoom,$
    widget,$
    pname0

  zoom = 1.0
  ;;-----------------------------------------------

  ;;----------------------------------------------------------
  ;; Return conditions:
  ;;
  if exists_kxkyspec eq 0 then return    
  if xregistered('spectrum_n_ave_see') then return
  ;;----------------------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)


  w1 = widget_button(row1,$
                     value='Plot',$
                     uvalue=1)

  w2 = widget_button(row1, $
                     value='ZOOM in', $
                     uvalue=2)

  w3 = widget_button(row1, $
                     value='ZOOM out', $
                     uvalue=3)

  w4 = widget_button(row1, $
                     value='PS dump', $
                     uvalue=4)

  w5 = widget_button(row1, $
                     value='data', $
                     uvalue=5)

  w6 = widget_button(row1,         $
                     value='Done', $
                     uvalue=6)

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

  xmanager,'spectrum_n_ave_see', $
    base, $
    event='spectrum_n_ave_event',$
    group_leader=group

end
