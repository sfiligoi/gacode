pro diffusion_i_ave_plot

  common GLOBAL
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA
  common ZMOMENT_DATA
  common PRIVATE_DIFFUSION_I_AVE

  ;;-------------------------------------------------------
  ;; Text labels
  ;;
  smf_tag,title,pname,ytitle

  title = title+' Diffusion'
  pname = 'diff_i-'+pname

  plot_def_new,pname

  print,format='(a,i2,a,i2,a)',$
     '(i_moment=',i_moment,', i_spec=',i_spec,')'

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
  y_i = fltarr(n_r,n_time)
  y_i_neo = fltarr(n_r,n_time)

  if i_f ge 0 then begin
     y_i[*,*] = diff_i[i_spec,i_f,i_moment,*,*]
  endif else begin
     y_i[*,*] = diff_i[i_spec,0,i_moment,*,*]+$
       diff_i[i_spec,1,i_moment,*,*]
  endelse

  y_i_neo[*,*] = 0.0
  if (neo_flag*exists_diff_i_ch eq 1) then begin
     if (i_spec eq 0 and i_moment eq 1) then begin
        y_i[*,*]     = y_i[*,*]+diff_i_ch[*,*]
        y_i_neo[*,*] = diff_i_ch[*,*]
     endif
  endif

  y       = fltarr(n_time)
  y_r     = fltarr(n_r)
  y_r_neo = fltarr(n_r)
  y_exp   = fltarr(n_r)

  for i=0,n_r-1 do begin
     ;; Set units here
     y[*] = y_i[i,*]*plot_units
     diff_stat_fast,y,it1,it2,ave_y
     y_r[i] = ave_y

     y[*] = y_i_neo[i,*]*plot_units
     diff_stat_fast,y,it1,it2,ave_y
     y_r_neo[i] = ave_y
     
  endfor
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; If the flag (i_div) is set, compute the conversion 
  ;; factor, gradient_factor, to convert the n' and T' in 
  ;; the definitions of chi/D to the full perturbed forms.
  ;;
  if (i_div eq 1) then begin


     gradient_factor = fltarr(n_r)

     y  = fltarr(n_r,n_time)
     y0 = fltarr(n_r)
     y1 = fltarr(n_r,n_time)

     if (i_moment eq 0) then begin

        ytitle = '(full dn/dr) '+ytitle

        y0[*] = dlnndr_s[i_spec,*]
        y1[*,*] = zmoment[*,i_spec,0,*]
        vector_deriv,r,y1,y,boundary_method
        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+den_s[i_spec,i]*y0[i]
           y[i,*] = y[i,*]/den_s[i_spec,i]
        endfor
     endif else begin

        ytitle = '(full dT/dr) '+ytitle

        y0[*] = dlntdr_s[i_spec,*]
        for i=0,n_r-1 do begin
           y1[i,*] = (zmoment[i,i_spec,1,*]- $
                      1.5*zmoment[i,i_spec,0,*]*tem_s[i_spec,i])/ $
             (1.5*(den_s[i_spec,i]+zmoment[i,i_spec,0,*]))        
        endfor
        vector_deriv,r,y1,y,boundary_method
        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+tem_s[i_spec,i]*y0[i]
           y[i,*] = y[i,*]/tem_s[i_spec,i]
        endfor
     endelse

     yt     = fltarr(n_time)
     y0_ave = fltarr(n_r)

     for i=0,n_r-1 do begin
        yt[*] = y[i,*]
        diff_stat_fast,yt,it1,it2,f_ave
        y0_ave[i] = f_ave
     endfor
     
     gradient_factor[*] = y0[*]/y0_ave[*]

     y_r[*] = y_r[*]*gradient_factor[*]

  endif 
  ;;-------------------------------------------------------
  
  ;;-------------------------------------------------------
  ;; Set the appropriate experimental chi-profile
  ;;
  y_exp[*] = 0.0

  if (i_moment eq 0) then begin
     if (i_spec eq n_kinetic-1) then begin
        y_exp[*] = diff_ne_exp[*]*plot_units
     endif
  endif else begin
     if (i_spec eq n_ion) then begin
        y_exp[*] = chi_e_exp[*]*plot_units
     endif else begin
        y_exp[*] = chi_i_exp[*]*plot_units
     endelse
  endelse

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
  if (neo_flag eq 1) then begin
   oplot,r,y_r_neo,color=color_vec[3]
  endif
  oplot,r,y_exp,color=color_vec[1]

  ;; Plot singular surface(s)

  get_singsurf_vec,n_ss,r_surf  
  n_surf = n_elements(r_surf)

  for i=0,n_surf-1 do begin
     oplot,r_surf[i]*[1,1],[-100,100],color=color_vec[3]
  endfor

  print,n_surf,r_surf

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
