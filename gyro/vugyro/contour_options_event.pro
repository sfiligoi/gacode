pro contour_options_event, contour_options

  common GLOBAL
  
  widget_control, contour_options.id, $
      get_uvalue=uvalue


  case (uvalue) of 
    
    1: begin
      c_table_min = contour_options.value
      print,'Set c_table_min =',c_table_min
    end

    2: begin
      c_table_max = contour_options.value      
      print,'Set c_table_max =',c_table_max
    end

    3: widget_control, contour_options.top, /destroy
    
  endcase

  return

end
