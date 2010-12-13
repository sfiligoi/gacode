pro transport_ne_te_ti_flows_see, group=group
;
; CH 6-19-08- minor bug fixes, added time-avg plot option
;

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_TR_FLOWS,$
    i_transport_ne_te_ti,$  
    i_ratio_value,$
    i_flow_avg_flag,$
    widget,$
    zoom

  zoom = 5.0
  i_ratio_value = 0
  i_flow_avg_flag = 0
  ;;-----------------------------------------------

  ;;-----------------------------------------------
  ;; Return conditions
  ;;
  if exists_transport_ne_te_ti eq 0 then return    
  if xregistered('transport_ne_te_ti_flows_see') then return
  ;;-----------------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  butt = widget_base(base,$
                     /row,$
                     /frame)

  w1 = widget_button(butt,$
                     value='Function',$
                     /menu)

  wlist = widget_button(w1,value='particle MW/Kev flows',uvalue=1)  
  wlist = widget_button(w1,value='Electron power flows',uvalue=2)  
  wlist = widget_button(w1,value='Ion power flows',uvalue=3)


  w4 = widget_button(butt,         $
                     value='+t', $
                     uvalue=4)

  w5 = widget_button(butt,         $
                     value='-t', $
                     uvalue=5)

  w6 = widget_button(butt,         $
                     value='++t', $
                     uvalue=6)

  w7 = widget_button(butt,         $
                     value='--t', $
                     uvalue=7)
  
  w8 = widget_button(butt,         $
                     value='all+++t', $
                     uvalue=8)

  w9 = widget_button(butt,         $
                     value='t-avg', $
                     uvalue=9)


  w14 = widget_button(butt, $
                     value='ZOOM in', $
                     uvalue=14)

  w15 = widget_button(butt, $
                     value='ZOOM out', $
                     uvalue=15)

  w12 =  widget_button(butt, $
                      value='ratio/value', $
                      uvalue=12)

  w11 = widget_button(butt, $
                      value='PS dump', $
                      uvalue=11)

  w13 = widget_button(butt, $
                      value='Done', $
                      uvalue=13)

  ;;---------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
      set_uvalue=state,$
      /no_copy, $
      /realize

  widget_control, draw, $
      get_value=widget

  xmanager,'transport_ne_te_ti_flows_see', $
      base, $
      event='transport_ne_te_ti_flows_event',$
      group_leader=group

end
