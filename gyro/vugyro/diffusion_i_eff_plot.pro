pro diffusion_i_eff_plot

  common GLOBAL
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA
  common ZMOMENT_DATA
  common PRIVATE_DIFFUSION_I_EFF

  ;;-------------------------------------------------------
  ;; Text labels
  ;;
  if (i_f eq 0) then begin
     tag_f = 'ES'
     ptag_f = 'es'
  endif else if (i_f eq 1) then begin
     tag_f = 'EM'
     ptag_f = 'em'
  endif else begin
     tag_f = 'Total'
     ptag_f = 'tot'
  endelse

  title = '!3Effective '+tag_f+' !4v!n'
  pname = 'diff_i_'+ptag_f+'_eff'
  ytitle = '!4v!3!deff!n = (!4v!3!di!n+!4v!3!de!n)/2'

  plot_def_new,pname

  if (i_units eq 0) then begin
     plot_units = 1.0
     units = ' [units of (c!ds!n/a)!4q!3!ds!n!u2!n]'
  endif else begin
     plot_units = xunits[8]
     units = ' [units of (m!u2!n/sec)]'
  endelse

  xtitle = '!3r/a with'+t_string
  title = title+' ['+strtrim(string(n_ss),2)+']'
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Pick out correct diffusivity and get the time 
  ;; average.
  ;;
  y_ion  = fltarr(n_field,n_r,n_time)
  y_elec = fltarr(n_field,n_r,n_time)

  y_ion[*,*,*] = 0.0
  for i=0,n_ion-1 do begin
     for j=0,n_time-1 do begin
        for k=0,n_field-1 do begin
           y_ion[k,*,j] = y_ion[k,*,j]+$
             den_s[i,*]*diff_i[i,k,1,*,j]/den_s[n_kinetic-1,*]
        endfor
     endfor
  endfor
  y_elec[*,*,*] = diff_i[n_kinetic-1,*,1,*,*]

  y_i = fltarr(n_r,n_time)
  if i_f ge 0 then begin
     y_i[*,*] = (y_ion[i_f,*,*]+y_elec[i_f,*,*])/2.0
  endif else begin
     y_i[*,*] = (y_ion[0,*,*]+y_elec[0,*,*])/2.0 + $
       (y_ion[1,*,*]+y_elec[1,*,*])/2.0
  endelse

  y     = fltarr(n_time)
  y_r   = fltarr(n_r)
  y_exp = fltarr(n_r)

  for i=1,n_r-1 do begin
     ;; Set units here
     y[*] = y_i[i,*]*plot_units
     diff_stat_fast,y,it1,it2,ave_y
     y_r[i] = ave_y
  endfor
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Set the appropriate experimental chi-profile
  ;;
  y_exp[*] = 0.5*(chi_i_exp[*]+chi_e_exp[*])*plot_units
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; PLOT results
  ;;
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

  oplot,r,y_r,color=color_vec[0]
  oplot,r,y_exp,color=color_vec[1]

  ;; Plot singular surface(s)

  get_singsurf_vec,n_ss,r_surf  
  n_surf = n_elements(r_surf)

  for i=0,n_surf-1 do begin
     oplot,r_surf[i]*[1,1],[-100,100],color=color_vec[3]
  endfor

  plot_bnd
  ;;-------------------------------------------------------

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_r-1 do begin
        printf,1,r[i],y_r[i]
     endfor
     printf,1,'Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
