pro exp_profiles_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_EXP_PROFILES,$
    i_profile,$
    i_r_rho,$
    widget

  i_r_rho = 0
  ;;-----------------------------------------------

  ;;------------------------------------------
  ;; Return conditions
  ;;
  if (exists_exp_profile eq 0) then return
  if xregistered('exp_profiles_see') then return
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

  for i=0,n_exp_profile-1 do begin
     x1 = widget_button(x,$
                        value=exp_profile_label[i],$
                        uvalue=i)  
  endfor

  x = widget_button(row1, $
                    value='-f/(df/dr)', $
                    uvalue=40)

  x = widget_button(row1, $
                    value='r/rho', $
                    uvalue=41)

  x = widget_button(row1, $
                    value='data', $
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

  xmanager,'exp_profiles_see', $
    base, $
    event='exp_profiles_event',$
    group_leader=group

end
