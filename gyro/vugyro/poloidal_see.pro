pro poloidal_see, group=group

  common GLOBAL
  common POLOIDAL_DATA

  if xregistered('poloidal_see') then return
  if exists_u eq 0 then return
  if simdir eq 'null' then return

  base = widget_base(title=simdir,$
                     /column)

  butt = widget_base(base,$
                     /row, $
                     /frame)

  w0 = widget_button(butt,$
                     value='<U>(t)',$
                     uvalue=0)

  w1 = widget_button(butt,$
                     value='U(t)',$
                     uvalue=1)

  w2 = widget_button(butt,$
                     value='U(n,t)',$
                     uvalue=15)

  w15 = widget_button(butt,$
                      value='u(n,t)',$
                      uvalue=2)

  w3 = widget_button(butt,         $
                     value='+t', $
                     uvalue=3)

  w4 = widget_button(butt,         $
                     value='-t', $
                     uvalue=4)

  w5 = widget_button(butt,         $
                     value='++t', $
                     uvalue=5)

  w6 = widget_button(butt,         $
                     value='--t', $
                     uvalue=6)

  w7 = widget_button(butt,         $
                     value='+n', $
                     uvalue=7)

  w8 = widget_button(butt,         $
                     value='-n', $
                     uvalue=8)

  w9 = widget_button(butt,         $
                     value='Pgrid', $
                     /menu)

  for i=0,n_mUt-1 do begin
    wlist = widget_button(w9,$
                          value=s_mUt[i],$
                          uvalue=20+i)  
  endfor

  w17 = widget_button(butt,         $
                      value='f tog',$
                      uvalue=17)

  ;;-----------------------------------

  current_val = widget_base(base,$
                            /row,$
                            /frame)

  w12 = widget_button(current_val, $
                      value='mpeg', $
                      uvalue=9)
  
  w13 = widget_button(current_val, $
                      value='PS dump', $
                      uvalue=12)
  
  w14 = widget_button(current_val, $
                      value='Done', $
                      uvalue=11)

  t_str = 't          '
  t_label = widget_label(current_val,value=t_str)

  n_str = 'n          '
  n_label = widget_label(current_val,value=n_str)

  jx_str = 'jx        '
  jx_label = widget_label(current_val,value=jx_str)

  ix_str = 'ix        '
  ix_label = widget_label(current_val,value=ix_str)

  state = {t_label:t_label, $
           n_label:n_label, $
           jx_label:jx_label, $
           ix_label:ix_label}

  ;;-----------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy*1.15)

  widget_control, base, $
      set_uvalue=state,$
      /no_copy, $
      /realize

  widget_control, draw, $
      get_value=poloidal_wid

  xmanager,'poloidal_see', $
      base, $
      event='poloidal_event',$
      group_leader=group

end
