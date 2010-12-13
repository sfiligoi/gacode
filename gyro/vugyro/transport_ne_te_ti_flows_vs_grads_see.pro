pro transport_ne_te_ti_flows_vs_grads_see, group=group

  common GLOBAL

  ;;-----------------------------------------------
  ;; Private (local) data  
  ;;
  common PRIVATE_TR_FLOWS_VS_GRADS,$
    i_transport_ne_te_ti,$  
    j_transport_ne_te_ti,$ 
    i_rho_offset,$
    widget,$
    zoom

  i_rho_offset = 0
  zoom = 5.0
  ;;-----------------------------------------------

  ;;-----------------------------------------------
  ;; Return conditions
  ;;
  if exists_transport_ne_te_ti eq 0 then return    
  if xregistered('transport_ne_te_ti_flows_vs_grads_see') then return
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



  w4 = widget_button(butt,$
                     value='-offset',$
                     uvalue=4)

  w5 = widget_button(butt,$
                     value='+off_set',$
                     uvalue=5)


  w14 = widget_button(butt, $
                      value='ZOOM in', $
                      uvalue=14)

  w15 = widget_button(butt, $
                      value='ZOOM out', $
                      uvalue=15)

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

  xmanager,'transport_ne_te_ti_flows_vs_grads_see', $
    base, $
    event='transport_ne_te_ti_flows_vs_grads_event',$
    group_leader=group

end
