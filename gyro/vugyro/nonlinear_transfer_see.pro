pro nonlinear_transfer_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data
  ;;
  common PRIVATE_NONLINEAR_TRANSFER,$
    TrEngGam_flag,$
    zoom_TrEngGam,$
    x2y_flag,$
    widget

  TrEngGam_flag = 0
  zoom_TrEngGam = 1.0
  x2y_flag = 0
  ;;-----------------------------------------------

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
  if n_n eq 1 then return
  if xregistered('nonlinear_transfer_see') then return
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
                    value='README', $
                    uvalue=11)

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=2)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=3)

  x = widget_button(row1, $
                    value='NLtrEnGam', $
                    uvalue=4)

  x = widget_button(row1, $
                    value='x2y shear', $
                    uvalue=5)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=6)

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
                  get_value=nonlinear_transfer_wid

  xmanager,'nonlinear_transfer_see', $
           base, $
           event='nonlinear_transfer_event',$
           group_leader=group

end
