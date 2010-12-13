pro neodiffusion_i_ave_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_NEODIFFUSION_I_AVE,$
    i_turb_neo,$
    i_div,$
    zoom,$
    widget

  i_turb_neo = 0
  i_div      = 0
  zoom       = 5.0
  ;;-----------------------------------------------

  ;;-----------------------------------------------
  ;; Return conditions
  ;;
  if (exists_diff_i eq 0) and (exists_diff_neo_i eq 0) then return 
  if xregistered('neodiffusion_i_ave_see') then return
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

  x = widget_button(row1,$
                    value='COMBINATION',$
                    /menu)

  tlevels=['Turbulent',$
           'Turb i+e energy',$
           'Turb+neo ion energy',$
           'Turb+neo electron energy',$
           'Turb+neo i+e energy',$
           'Neo only']

  for i=0,5 do begin
     x1 = widget_button(x,$
                        value=tlevels[i],$
                        uvalue=30+i)  
  endfor
  
  x = widget_button(row1,$
                    value='DIVISOR',$
                    /menu)

  tlevels=['NORMAL','LOCAL DIFF ','DIV FLOW']
  for i=0,2 do begin
     x1 = widget_button(x,$
                        value=tlevels[i],$
                        uvalue=20+i)  
  endfor

  x = widget_button(row1, $
                    value='ZOOM in', $
                    uvalue=9)

  x = widget_button(row1, $
                    value='ZOOM out', $
                    uvalue=10)

  x = widget_button(row1, $
                    value='zero', $
                    uvalue=5)

  x =  widget_button(row1, $
                     value='units', $
                     uvalue=11)

  row2 = widget_base(base, $
                     /row, $
                     /frame)

  x = widget_button(row2, $
                    value='ss+', $
                    uvalue=12)

  x = widget_button(row2, $
                    value='ss-', $
                    uvalue=13)

  x = widget_button(row2, $
                    value='data', $
                    uvalue=14)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=15)

  x = widget_button(row2, $
                    value='Done', $
                    uvalue=16)

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

  xmanager,'neodiffusion_i_ave_see', $
    base, $
    event='neodiffusion_i_ave_event',$
    group_leader=group

end
