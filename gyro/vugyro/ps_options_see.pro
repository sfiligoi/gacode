pro ps_options_see, group=group

   common GLOBAL

   if xregistered('ps_options_see') then return

   base = widget_base(title='Postscript Options',$
                      /column)

   p1 = cw_field(base,$
                 /frame,$
                 /row,$
                 title='xsize: ',$
                 uvalue=1,$
                 xsize=5, $
                 value=xsize,$
                 /return_events)

   p2 = cw_field(base,$
                 /frame,$
                 /row,$
                 title='ysize: ',$
                 uvalue=2,$
                 xsize=5, $
                 value=ysize,$
                 /return_events)

   p3 = cw_bgroup(base,$
                  ['ps','eps'],$
                  /row, $
                  /exclusive, $
                  label_top='Format',$
                  /frame, $
                  uvalue=3, $
                  set_value=ps_val)

   p3 = cw_bgroup(base,$
                  ['b/w','color'],$
                  /row, $
                  /exclusive, $
                  label_top='Drawing',$
                  /frame, $
                  uvalue=4, $
                  set_value=ps_color)

   p4 = cw_bgroup(base,$
                  ['on','off'],$
                  /row, $
                  /exclusive, $
                  label_top='Sim label',$
                  /frame, $
                  uvalue=5, $
                  set_value=ps_label)

   p5 = cw_bgroup(base,$
                  ['on','off'],$
                  /row, $
                  /exclusive, $
                  label_top='Thin lines',$
                  /frame, $
                  uvalue=6, $
                  set_value=ps_thick)

   p6 = widget_button(base,         $
                      value='Done', $
 		      uvalue=7)

   widget_control, base, $
                   /no_copy, $
                   /realize

   xmanager,'ps_options_see', $
            base, $
	    event='ps_options_event',$
            group_leader=group

end
