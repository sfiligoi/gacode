pro master_see, group=group

  common GLOBAL

  ;;---------------------------------------------------------
  ;; Return conditions
  ;;
  if simdir eq 'null' then return
  if xregistered('master_see') then return
  ;;----------------------------------------------------------

  base = widget_base(title='MASTER',$
                     /row)

  ;;----------------------------------------------------------
  ;; BUTTONS
  ;;----------------------------------------------------------

  ;; Column 1

  col1 = widget_base(base,$
                     /column, $
                     /frame)

  x = cw_field(col1,$
               /return_events, $
               title='t_min', $
               /floating, $
               value=0.0, $
               xsize=8, $
               uvalue=1)

  x = cw_field(col1,$
               /return_events, $
               title='t_max', $
               /floating, $
               value=0.0, $
               xsize=8, $
               uvalue=2)

  x = widget_button(col1, $
                    value='Color Table', $
                    /menu)

  for i=0,n_CTab-1 do begin
     x1 = widget_button(x,$
                        value=s_CTab[i],$
                        uvalue=50+i)  
  endfor
  
  x = widget_button(col1, $
                    value='Levels', $
                    /menu)

  for i=0,n_nLevels-1 do begin
     x1 = widget_button(x,$
                        value=s_nLevels[i],$
                        uvalue=40+i)  
  endfor

  pair = widget_base(col1,$
                     /row, $
                     /frame)

  x1 = widget_button(pair, $
                     value='t=tmin', $
                     uvalue=3)

  x1 = widget_button(pair, $
                     value='t=tmax', $
                     uvalue=4)

  pair = widget_base(col1,$
                     /row, $
                     /frame)

  x1 = widget_button(pair, $
                     value='bin+', $
                     uvalue=5)

  x1 = widget_button(pair, $
                     value='bin-', $
                     uvalue=6)

  x = widget_button(col1,         $
                    value='Done', $
                    uvalue=7)

  ;; Column 2

  col2 = widget_base(base,$
                     /column, $
                     /frame)

  t_label = widget_label(col2,$
                         value='Average Time:'+t_string+' ',$
                         /align_left)

  ct_label = widget_label(col2,$
                          value=' Color Table: '+s_CTab[CTab]+'          ',$
                          /align_left)

  level_str = strtrim(string(nLevels),2)

  cl_label = widget_label(col2,$
                          value='      Levels: '+level_str+'   ',$
                          /align_left)

  cnt_str = strtrim(string(t_c),2)

  cnt_label = widget_label(col2,$
                           value='         t_c: '+cnt_str+'   ',$
                           /align_left)

  bin_str = strtrim(string(nbin),2)

  bin_label = widget_label(col2,$
                           value='        nbin: '+bin_str+'   ',$
                           /align_left)

  state = { t_label:t_label,$
            ct_label:ct_label,$
            cl_label:cl_label,$
            cnt_label:cnt_label, $
            bin_label:bin_label }

  widget_control, base, $
                  set_uvalue=state,$
                  /no_copy, $
                  /realize

  xmanager,'master_see', $
           base, $
           event='master_event',$
           group_leader=group

end
