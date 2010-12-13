pro contour_options_see, group=group

   common GLOBAL

   if xregistered('contour_options_see') then return
   if simdir eq 'null' then return

   base = widget_base(title=simdir,$
                      /column)

   p1 = cw_field(base,$
                 /frame,$
                 /row,$
                 /floating,$
                 title='min: ',$
                 xsize=10, $
                 value=c_table_min,$
                 uvalue=1,$
                 /return_events)

   p2 = cw_field(base,$
                 /frame,$
                 /row,$
                 /floating,$
                 title='max: ',$
                 xsize=10, $
                 value=c_table_max,$
                 uvalue=2,$
                 /return_events)

   p3 = widget_button(base,         $
                      value='Done', $
 		      uvalue=3)

   widget_control, base, $
                   /no_copy, $
                   /realize

   xmanager,'contour_options_see', $
            base, $
	    event='contour_options_event',$
            group_leader=group

end
