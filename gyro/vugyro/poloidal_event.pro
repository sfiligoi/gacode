pro poloidal_event, poloidal

  common GLOBAL
  common PLOT_VARIABLES
  common POLOIDAL_DATA

  widget_control, poloidal.id, $
      get_uvalue=uvalue

  wset, poloidal_wid
  
  ;; poloidal refinement
  if (uvalue ge 20) then begin
    mUt = v_mUt[uvalue-20]
    goto, plot_it
  endif

  top:

  ;;--------------------------
  ;; Now, begin menu selection
  ;;
  case (uvalue) of 
    
    0: begin   
      active_poloidal = 10
      ave_poloidal = 1

      set_contour_levels

      one_contour = 0
      goto, plot_it
    end

    1: begin
      active_poloidal = 10
      ave_poloidal = 0

      set_contour_levels

      one_contour = 0
      goto, plot_it
    end

    15: begin
      active_poloidal = 10
      ave_poloidal = 0
      set_contour_levels
      one_contour = 1
      goto, plot_it
    end

    2: begin

      active_poloidal = 2

      set_contour_levels
      ny = n_theta_plot*mUt

      plot_def_new,'poloidal_orig'+strtrim(string(t_c),2)

      ;;-----------------------------------
      ;; widget control for label functions
      ;;
      widget_control, poloidal.top, $
          get_uvalue=state, $
          /no_copy

      t_str = 't = '+strtrim(string(t[t_c]),2)
      n_str = 'n = '+strtrim(string(n_tor[in_c]),2)
      jx_str = 'jx = '+strtrim(string(fix(ny)),2)
      ix_str = 'ix = '+strtrim(string(fix(n_r)),2)

      widget_control, state.t_label, set_value=t_str
      widget_control, state.n_label, set_value=n_str      
      widget_control, state.jx_label, set_value=jx_str
      widget_control, state.ix_label, set_value=ix_str      

      widget_control, poloidal.top, $
          set_uvalue=state, $
          /no_copy
      ;;-------------------------------------

      a = fltarr(ny,n_r)
      U_temp = fltarr(n_theta_plot,n_r)

      U_temp[*,*] = U[0,*,*,0,in_c,t_c]
      a = rebin(U_temp,ny,n_r)

      make_fine_grid
      torus_contour,t_c

    end   

    
    3: begin
      counter_up,t_c,n_time1,1
      goto, plot_it
    end

    4: begin
      counter_dn,t_c,0,1
      goto, plot_it
    end

    5: begin
      counter_up,t_c,n_time1,10 
      goto, plot_it
    end

    6: begin 
      counter_dn,t_c,0,10
      goto, plot_it
    end

    7: begin
      counter_up,in_c,n_n-1,1 
      goto, plot_it
    end

    8: begin 
      counter_dn,in_c,0,1
      goto, plot_it
    end
    
    9: begin
      
      plot_def_new,'poloidal_orig0'

      make_fine_grid

      openw,1,'incl'
      
      line = 1

      for tt=0,n_time1 do begin
        
        if (t[tt] ge t_min) and (t[tt] le t_max) then begin

          if one_contour eq 0 then begin
            extract_contour,tt
          endif else begin
            extract_one_contour,in_c,tt
          endelse

          torus_contour_movie,tt

          file = 'torus'+strtrim(string(tt),2)+'.jpg'
          tmp_image = tvrd(true=3,/order) 
          write_jpeg,file,tmp_image,true=3,/order
          print,'Wrote '+file
          printf,1,file

        endif

      endfor
      
      close,1

    end

    10: begin 

      active_poloidal=10

      ny = n_theta_plot*mUt

      ;; fix scaling
      pname='contour'
      plot_def_new,pname

      ;;-------------------------------------------------
      ;; widget control for label functions
      ;;
      widget_control, poloidal.top, $
          get_uvalue=state, $
          /no_copy

      t_str = 't = '+strtrim(string(t[t_c]),2)
      if one_contour eq 0 then begin
        n_str = 'n = (all)'
      endif else begin
        n_str = 'n = '+strtrim(string(n_tor[in_c]),2)
      endelse
      jx_str = 'jx = '+strtrim(string(fix(ny)),2)
      ix_str = 'ix = '+strtrim(string(fix(n_r)),2)

      widget_control, state.t_label, set_value=t_str
      widget_control, state.n_label, set_value=n_str      
      widget_control, state.jx_label, set_value=jx_str
      widget_control, state.ix_label, set_value=ix_str      

      widget_control, poloidal.top, $
          set_uvalue=state, $
          /no_copy
      ;;-------------------------------------------------
      
      make_fine_grid

      if one_contour eq 0 then begin
        extract_contour,t_c
      endif else begin
        extract_one_contour,in_c,t_c
      endelse

      torus_contour,t_c

    end

    11: widget_control, poloidal.top, /destroy     

    12: begin
      plot_mode = 2
      goto, plot_it
    end

    17: begin
      if (n_field gt 1) then i_f = 1-i_f
      goto, plot_it
    end

  endcase
  ;;
  ;;----------------------------------------------------------

  return

  plot_it:

  ;; define top label

  if (one_contour eq 0) then begin

    if (i_f eq 0 and ave_poloidal eq 1) then begin
      title = "!4u!3-<!4u!3>" 
    endif

    if (i_f eq 0 and ave_poloidal eq 0) then begin
      title = "!4u!3" 
    endif

    if (i_f eq 1 and ave_poloidal eq 1) then begin
      title = "!3A!d!9#!3!n-<!3A!d!9#!3!n>" 
    endif

    if (i_f eq 1 and ave_poloidal eq 0) then begin
      title = "!3A!d!9#!3!n"
    endif

  endif else begin

    if (i_f eq 0) then begin
      title = "!4u!3!dn!n" 
    endif

    if (i_f eq 1) then begin
      title = "!3A!d!9#!3!n!dn!n" 
    endif
  endelse
 
  uvalue=active_poloidal
  goto,top

end
