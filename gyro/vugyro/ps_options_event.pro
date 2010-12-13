pro ps_options_event, ps_options

   common GLOBAL
   
   widget_control, ps_options.id, $
                   get_uvalue=uvalue

   case (uvalue) of 
   
      1: xsize = ps_options.value
      
      2: ysize = ps_options.value    

      3: ps_val = ps_options.value

      4: ps_color = ps_options.value

      5: ps_label = ps_options.value

      6: ps_thick = ps_options.value

      7: widget_control, ps_options.top, /destroy
      
   endcase

   return

end
