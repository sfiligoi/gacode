pro master_event, master

  common GLOBAL
  common PLOT_VARIABLES
  
  widget_control, master.id, $
                  get_uvalue=uvalue

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  ;; Color tables
  if (uvalue ge 50) then begin     
     CTab = v_CTab[uvalue-50]
     ;;device,true_color=24,decomposed=0
     loadct,CTab
     goto, plot_it
  endif

  ;; Color contour levels
  if (uvalue ge 40) then begin
     nLevels = v_nLevels[uvalue-40]
     goto, plot_it
  endif

  case (uvalue) of 
     
     1: begin
        t_min = master.value
        goto, plot_it
     end

     2: begin
        t_max = master.value
        goto, plot_it
     end

     3: begin
        t_c = 0
        goto, plot_it
     end

     4: begin
        t_c = n_time1
        goto, plot_it
     end

     5: begin
        nbin = nbin+1
        goto, plot_it
     end

     6: begin
        nbin = nbin-1
        if (nbin lt 2) then nbin=2
        goto, plot_it
     end

     7: widget_control, master.top, /destroy
     
  endcase

  return

  plot_it:

  widget_control, master.top, $
                  get_uvalue=state, $
                  /no_copy

  ;; reset time-interval string 

  get_t_string
  t_indices,t_min,t_max,it1,it2

  widget_control, $
    state.t_label, $
    set_value='Average Time:'+t_string

  widget_control, $
    state.ct_label, $
    set_value=' Color Table: '+s_CTab[CTab]

  level_str = strtrim(string(nLevels),2)
 
  widget_control, $
    state.cl_label, $
    set_value='      Levels: '+level_str

  cnt_str   = strtrim(string(t_c),2)

  widget_control, $
    state.cnt_label, $
    set_value='         t_c: '+cnt_str

  bin_str = strtrim(string(nbin),2)

  widget_control, $
    state.bin_label, $
    set_value='        nbin: '+bin_str

  widget_control, master.top, $
                  set_uvalue=state, $
                  /no_copy

  return

end
