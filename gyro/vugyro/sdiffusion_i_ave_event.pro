pro sdiffusion_i_ave_event, sdiffusion_i_ave

  common GLOBAL
  common TAG_DATA
  common PLOT_VARIABLES
  common PRIVATE_SDIFFUSION_I_AVE
  
  widget_control, sdiffusion_i_ave.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: begin
        i_spec = i_spec+1
        if (i_spec ge n_kinetic) then i_spec = 0
        goto, plot_it
     end   

     2: begin
        i_moment = 1+i_moment
        if (i_moment ge 2) then i_moment = 0
        goto, plot_it
     end   

     3: begin
        zoom = zoom/2
        goto, plot_it
     end

     4: begin
        zoom = zoom*2
        goto, plot_it
     end

     5: begin
        i_units = 1-i_units
        goto, plot_it
     end

     6: begin
        if (exists_sdiff_i_t eq 1) then begin
           i_tp = i_tp+1
           if i_tp ge 3 then i_tp = 0
        endif
        goto, plot_it
     end

     7: begin
        i_zero = 1-i_zero
        goto, plot_it
     end

     8: begin
        plot_mode = 2
        goto, plot_it
     end

     9: widget_control, sdiffusion_i_ave.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  y_i = fltarr(n_r,n_time)

  y_i_1 = fltarr(n_r,n_time)
  y_i_2 = fltarr(n_r,n_time)
  y_i_3 = fltarr(n_r,n_time)
  y_i_4 = fltarr(n_r,n_time)

  
  tag_gen,n_ion,i_spec,i_moment,0,0

  title = tag_s
  pname = ftag_s
  if (i_moment eq 0) then begin
     title = title+' Mom.'
     pname = pname+'_mom'
  endif else begin
     title = title+' Heating'
     pname = pname+'_heat'
  endelse
  ytitle = pname

  case (i_tp) of

     0: begin

        pname = 'sdiff_i-'+pname

        y_i[*,*] = sdiff_i[i_spec,i_moment,*,*]

       y_i_1[*,*] = spdiff_i[i_spec,i_moment,0,*,*]
       y_i_2[*,*] = spdiff_i[i_spec,i_moment,1,*,*]
       y_i_3[*,*] = spdiff_i[i_spec,i_moment,2,*,*]

       ;;  y_i_4[*,*] = y_i_1[*,*]+y_i_2[*,*]+y_i_3[*,*] 

       ;; y_i[*,*] = y_i_4[*,*]




     end

     1: begin

        pname = 'sdiff_i_t-'+pname
        title = title+' (trap)'
        y_i[*,*] = sdiff_i_trapped[i_spec,i_moment,*,*]

     end

     2: begin

        pname = 'sdiff_i_p-'+pname
        title = title+' (pass)'
        y_i[*,*] = sdiff_i[i_spec,i_moment,*,*]-$
          sdiff_i_trapped[i_spec,i_moment,*,*]

     end

  endcase

  if (i_units eq 0) then begin
     plot_units = 1.0
     units = ' [units of '+chi_string+']'
  endif else begin
     plot_units = xunits[8]
     units = ' [units of (m!u2!n/sec)]'
  endelse

  y_exp = fltarr(n_r)
  y2_exp = fltarr(n_r)
  y_exp[*]=0.0
  y2_exp[*]=0.0

  if (i_moment eq 0) then begin
     y_exp[*] = eta_i_tot_exp[*]*plot_units
     y2_exp[*] = eta_i_diff_exp[*]*plot_units
  endif

  ;; Diffusion TIME trace

  plot_def_new,pname

  y       = fltarr(n_time)
  A_f_i   = fltarr(n_r)

  for ii=n_bnd,n_r-1-n_bnd do begin
     ;; Set units here
     y[*] = y_i[ii,*]*plot_units
     diff_stat_fast,y,it1,it2,A_f
     A_f_i[ii] = A_f
  endfor

  if i_tp eq 0 then begin

  y1       = fltarr(n_time)
  A_f_i_1   = fltarr(n_r)

  for ii=n_bnd,n_r-1-n_bnd do begin
     ;; Set units here
     y1[*] = y_i_1[ii,*]*plot_units
     diff_stat_fast,y1,it1,it2,A_f
     A_f_i_1[ii] = A_f
  endfor

  y2      = fltarr(n_time)
  A_f_i_2   = fltarr(n_r)

  for ii=n_bnd,n_r-1-n_bnd do begin
     ;; Set units here
     y2[*] = y_i_2[ii,*]*plot_units
     diff_stat_fast,y2,it1,it2,A_f
     A_f_i_2[ii] = A_f
  endfor


  y3      = fltarr(n_time)
  A_f_i_3   = fltarr(n_r)

  for ii=n_bnd,n_r-1-n_bnd do begin
     ;; Set units here
     y3[*] = y_i_3[ii,*]*plot_units
     diff_stat_fast,y3,it1,it2,A_f
     A_f_i_3[ii] = A_f
  endfor

  endif

  xtitle = '!3r/a with'+t_string

  plot,[0],[0],$
    /nodata,$
    title=title,$
    xstyle=1,$
    xrange=[min(r),max(r)],$
    xtitle=xtitle,$
    ystyle=1,$
    yrange=[-i_zero*zoom,zoom],$
    ytitle=ytitle+units,$
    color=line

  oplot,r,A_f_i,color=color_vec[0]
    ;; 0=blue  total
  if i_tp eq 0 then begin
   oplot,r,A_f_i_1,color=color_vec[1]
    ;; 1=purple M=z H=par
   oplot,r,A_f_i_2,color=color_vec[2]
    ;; 2=yellow M=y H=perp
   oplot,r,A_f_i_3,color=color_vec[3]
    ;; 3=green  M=u H=0
  endif 
  oplot,r,y_exp,color=color_vec[1]
  oplot,r,y2_exp,color=color_vec[2]

  plot_bnd
  plot_finish

  return

end
