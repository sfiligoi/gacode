pro collision_see, group=group

  common GLOBAL

  ;;-------------------------------------
  ;; Return conditions
  ;;
  if exists_coll eq 0 then return
  if xregistered('collision_see') then return
  ;;-------------------------------------

  base = widget_base(title=simdir,$
                     /column)

  row1 = widget_base(base,$
                     /row)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  ;; Row 1

  x = widget_button(row1,$
                    value='Plot',$
                    uvalue=1)

  x = widget_button(row1,         $
                    value='+m', $
                    uvalue=2)

  x = widget_button(row1,         $
                    value='-m', $
                    uvalue=3)

  x = widget_button(row1,         $
                    value='+k', $
                    uvalue=4)

  x = widget_button(row1,         $
                    value='-k', $
                    uvalue=5)

  x = widget_button(row1,         $
                    value='+i', $
                    uvalue=6)

  x = widget_button(row1,         $
                    value='-i', $
                    uvalue=7)

  ;; Row 2

  row2 = widget_base(base,$
                     /row,$
                     /frame)

  x = widget_button(row2, $
                    value='PS dump', $
                    uvalue=8)

  x = widget_button(row2,         $
                    value='Done', $
                    uvalue=9)

  m_str = ' (null)        '
  m_label = widget_label(row2,value=m_str)

  k_str = ' (null)        '
  k_label = widget_label(row2,value=k_str)

  i_str = ' (null)     '
  i_label = widget_label(row2,value=i_str)

  state = {m_label:m_label, $
           k_label:k_label, $
           i_label:i_label }

  ;;----------------------------------------------------------
  ;; DRAW WIDGET and CONTROL
  ;;----------------------------------------------------------

  draw = widget_draw(base,     $
                     xsize=sx, $
                     ysize=sy)

  widget_control, base, $
                  set_uvalue=state, $
                  /no_copy, $
                  /realize

  widget_control, draw, $
                  get_value=collision_wid

  xmanager,'collision_see', $
           base, $
           event='collision_event',$
           group_leader=group

end
