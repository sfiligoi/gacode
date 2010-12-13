pro all_profiles_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_ALL_PROFILES,$
    i_all_profile,$
    widget
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if simdir eq 'null' then return
  if xregistered('all_profiles_see') then return
  ;;------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  row1 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row1,$
                    value='SELECT PROFILE',$
                    /menu)

  for i=0,n_profile_label-1 do begin
     x1 = widget_button(x,$
                        value=profile_label[i],$
                        uvalue=i+1)  
  endfor

  x = widget_button(row1, $
                    value='i/e', $
                    uvalue=40)

  x = widget_button(row1, $
                    value='f/(df/dr)', $
                    uvalue=41)

  x = widget_button(row1, $
                    value='export', $
                    uvalue=42)

  x = widget_button(row1, $
                    value='PS dump', $
                    uvalue=43)

  x = widget_button(row1, $
                    value='Done', $
                    uvalue=44)

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
    /no_copy, $
    /realize

  widget_control, draw, $
    get_value=widget

  xmanager,'all_profiles_see', $
    base, $
    event='all_profiles_event',$
    group_leader=group

end
