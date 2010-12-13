pro flow_i_ave_event, flow_i_ave

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_FLOW_I_AVE
  
  widget_control, flow_i_ave.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------
  
  case (uvalue) of 
     
     1: begin
        active_flow_i_ave = 1
        title = 'Ion Energy Flow'
        pname = 'ion_energy_flow'
        ytitle = '!3Flow in MW'
        goto, plot_it
     end   

     2: begin
        active_flow_i_ave = 2
        title = 'Electron Energy Flow'
        pname = 'elec_energy_flow'
        ytitle = '!3Flow in MW'
        goto, plot_it
     end   

     3: begin
        active_flow_i_ave = 3 
        title = 'Electron Plasma Flow'
        pname = 'electron_plasma_flow'
        ytitle = '!3Flow in MW/Kev'
        goto, plot_it
     end

     4: begin
        active_flow_i_ave = 4 
        title = 'Angular Momentum Flow'
        pname = 'angular_momentum_flow'
        ytitle = '!3Flow in Newtons-meters'
        goto, plot_it
     end

     5: begin
        active_flow_i_ave = 5
        title = 'Ion Heating Flow'
        pname = 'ion_heating_flow'
        ytitle = '!3Flow in MW'
        goto, plot_it
     end

     6: begin
        active_flow_i_ave = 6
        title = 'Electron Heating Flow'
        pname = 'electron_heating_flow'
        ytitle = '!3Flow in MW'
        goto, plot_it
     end

     7: begin
        active_flow_i_ave = 7
        title = 'Total Heating Flow'
        pname = 'total_heating_flow'
        ytitle = '!3Flow in MW'
        goto, plot_it
     end

     8: begin
        zoom = zoom/2 
        goto, plot_it
     end

     9: begin
        zoom = zoom*2 
        goto, plot_it
     end

     91: begin
        i_zero = 1-i_zero 
        goto, plot_it
     end


     10: begin
        plot_mode = 2
        goto, plot_it
     end

     11: begin
        plot_export = 1
        goto, plot_it
     end

     111: begin
       neo_flag = 1-neo_flag
       goto, plot_it
     end

     12: widget_control, flow_i_ave.top, /destroy
     
  endcase
  
  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  plot_def_new,pname

  y_i   = fltarr(n_r,n_time)
  y_i_neo   = fltarr(n_r,n_time)
  y_i_add = fltarr(n_r,n_time)
  y_t   = fltarr(n_time)
  y     = fltarr(n_r)
  y_neo = fltarr(n_r)
  y_exp = fltarr(n_r)
  y2_exp = fltarr(n_r)

  y_exp[*]=0.0
  y2_exp[*]=0.0

  case (active_flow_i_ave) of

     1: begin
        for i=0,n_r-1 do begin
           y_i[i,*] = diff_i[0,0,1,i,*]*diff_to_flow_e1[i]
        endfor

        if (n_field eq 2) then begin
           for i=0,n_r-1 do begin
              y_i[i,*] = y_i[i,*]+diff_i[0,1,1,i,*]*diff_to_flow_e1[i]
           endfor
        endif

      y_i_neo[*,*] = 0.
      if (exists_diff_i_ch eq 1) then begin
        if (neo_flag eq 1) then begin
           for i=0,n_r-1 do begin
              y_i[i,*] = y_i[i,*]+diff_i_ch[i,*]*diff_to_flow_e1[i]
              y_i_neo[i,*] = diff_i_ch[i,*]*diff_to_flow_e1[i]
           endfor
        endif
      endif

        for i=n_bnd,n_r-1-n_bnd do begin
           y_exp[i] = chi_i_exp[i]*diff_to_flow_e1[i]
        endfor
     end

     2: begin
        for i=0,n_r-1 do begin
           y_i[i,*] = diff_i[n_kinetic-1,0,1,i,*]*diff_to_flow_e2[i]
        endfor

        if (n_field eq 2) then begin
           for i=0,n_r-1 do begin
              y_i[i,*] = y_i[i,*]+diff_i[n_kinetic-1,1,1,i,*]*diff_to_flow_e2[i]
           endfor
        endif

        for i=n_bnd,n_r-1-n_bnd do begin
           y_exp[i] = chi_e_exp[i]*diff_to_flow_e2[i]
        endfor
     end

     3: begin
        for i=0,n_r-1 do begin
           y_i[i,*] = diff_i[n_kinetic-1,0,0,i,*]*diff_to_flow_ne[i]
        endfor
        if (n_field eq 2) then begin
           for i=0,n_r-1 do begin
              y_i[i,*] = y_i[i,*]+diff_i[n_kinetic-1,1,0,i,*]*diff_to_flow_ne[i]
           endfor
        endif

        for i=n_bnd,n_r-1-n_bnd do begin
           y_exp[i] = diff_ne_exp[i]*diff_to_flow_ne[i]
        endfor
     end

     4: begin
        for i=0,n_r-1 do begin
           y_i[i,*] = sdiff_i[0,0,i,*]*diff_to_flow_mi[i]
        endfor
        for i=n_bnd,n_r-1-n_bnd do begin
           y_exp[i] = eta_i_tot_exp[i]*diff_to_flow_mi[i]
           y2_exp[i] = eta_i_diff_exp[i]*diff_to_flow_mi[i]
        endfor
     end


     5: begin
        y_i_add[*,*] = 0.0
        y_i[*,*] = 0.0
        for i=n_bnd,n_r-1-n_bnd do begin
           for is=0,n_kinetic-2 do begin
              y_i_add[i,*] =  y_i_add[i,*]+sdiff_i[is,1,i,*]*diff_to_flow_heating[i]*(r(i)-r(i-1))
           endfor
           y_i[i,*] = y_i[i-1,*] +y_i_add[i,*] 
        endfor
     end

     6: begin
        y_i[*,*] = 0.0
        for i=n_bnd,n_r-1-n_bnd do begin
           y_i[i,*] = y_i[i-1,*]+ sdiff_i[n_kinetic-1,1,i,*]*diff_to_flow_heating[i]*(r(i)-r(i-1))
        endfor
     end


     7: begin
        y_i_add[*,*] = 0.0
        y_i[*,*] = 0.0
        for i=n_bnd,n_r-1-n_bnd do begin
           for is=0,n_kinetic-1 do begin
              y_i_add[i,*] =  y_i_add[i,*]+sdiff_i[is,1,i,*]*diff_to_flow_heating[i]*(r(i)-r(i-1))
           endfor
           y_i[i,*] = y_i[i-1,*] +y_i_add[i,*]
        endfor

     end


  endcase

  for i=0,n_r-1 do begin
     y_t[*] = y_i[i,*]
     diff_stat_fast,y_t,it1,it2,ave_y
     y[i] = ave_y

     y_t[*] = y_i_neo[i,*]
     diff_stat_fast,y_t,it1,it2,ave_y
     y_neo[i] = ave_y

  endfor

  xtitle = '!3r/a with '+t_string

  plot,[0],[0],$
    /nodata,$
    title=title,$
    xstyle=1,$
    xrange=[min(r),max(r)],$
    xtitle=xtitle,$
    ystyle=1,$
    yrange=[-i_zero*zoom,zoom],$
    ytitle=ytitle,$
    color=line

  oplot,r,y,color=color_vec[0]
  if(neo_flag eq 1) then begin
    oplot,r,y_neo,color=color_vec[3]
  endif
  oplot,r,y_exp,color=color_vec[1]
  oplot,r,y2_exp,color=color_vec[2]

  plot_bnd 
  plot_finish

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_r-1 do begin
        printf,1,r[i],y[i],y_exp[i]
     endfor
     printf,1,'Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------
  
  return

end
