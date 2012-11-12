pro osc_plot
; CH 7.31.09- fixed incorrect indexing on line 75 of z_moment
;
  common GLOBAL
  common PLOT_VARIABLES
  common TAG_DATA
  common ZMOMENT_DATA
  common PROFILE_SIM_DATA
  common PRIVATE_OSC

  tag_gen,n_ion,i_spec,i_moment,i_f,1

  y    = fltarr(n_r,n_time)
  y_eq = fltarr(n_r)

  if i_gradient eq 0 then begin 

     if (i_moment eq 0) then begin
        title  = tag_s+' '+tag_m  
        ytitle = '<n'+ftag_s+'>'
        pname  = 'osc-'+ftag_m+'-'+ftag_s

        y_eq[*] = equil_flag*den_s[i_spec,*]

        for i=0,n_r-1 do begin
           y[i,*] = y_eq[i]+zmoment[i,i_spec,0,*]
        endfor
     endif else begin 
        title  = tag_s+' '+tag_m  
        ytitle = '<T'+ftag_s+'>'
        pname  = 'osc-'+ftag_m+'-'+ftag_s

        y_eq[*] = equil_flag*tem_s[i_spec,*]

        for i=0,n_r-1 do begin
           y[i,*] = y_eq[i]+(zmoment[i,i_spec,1,*]- $
                             1.5*zmoment[i,i_spec,0,*]*tem_s[i_spec,i])/ $
             (1.5*(den_s[i_spec,i]+zmoment[i,i_spec,0,*]))
        endfor
     endelse

  endif else begin

     if (i_div eq 0) then begin
        div_tag = ''
     endif else begin
        div_tag = ' (full T)'
     endelse

     if (i_moment eq 0) then begin
        title  = tag_s+' '+tag_m+' Gradient'  
        ytitle = '<dlnn'+ftag_s+'/dr>'+div_tag
        pname  ='osc-d'+ftag_m+'-'+ftag_s

        y_eq[*] = dlnndr_s[i_spec,*]

        yx = fltarr(n_r,n_time)
        yx[*,*] = zmoment[*,i_spec,0,*]
        
        vector_deriv,r,yx,y,boundary_method

        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+den_s[i_spec,i]*y_eq[i]
           y[i,*] = y[i,*]/den_s[i_spec,i]
        endfor
     endif else begin 
        title  = tag_s+' '+tag_m+' Gradient'  
        ytitle = '<dlnT'+ftag_s+'/dr>'+div_tag
        pname  = 'osc-d'+ftag_m+'-'+ftag_s

        y_eq[*] = dlntdr_s[i_spec,*]

        yx = fltarr(n_r,n_time)

        for i=0,n_r-1 do begin
           yx[i,*] = (zmoment[i,i_spec,1,*]- $
                      1.5*zmoment[i,i_spec,0,*]*tem_s[i_spec,i])/ $
             (1.5*(den_s[i_spec,i]+zmoment[i,i_spec,0,*]))        
        endfor

        vector_deriv,r,yx,y,boundary_method

        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+tem_s[i_spec,i]*y_eq[i]
           y[i,*] = y[i,*]/tem_s[i_spec,i]
        endfor
     endelse

  endelse

  plot_def_new,pname

  xtitle = '!3r/a with'+t_string

  y0     = fltarr(n_time)
  y0_ave = fltarr(n_r)

  for i=0,n_r-1 do begin
     y0[*] = y[i,*]
     diff_stat_fast,y0,it1,it2,f_ave
     y0_ave[i] = f_ave
  endfor
  
  y_mid = y_eq[n_r/2]
  y_min=min(y0_ave[n_bnd:n_r-1-n_bnd])
  y_max=max(y0_ave[n_bnd:n_r-1-n_bnd])
  y_mid=0.5*(y_min+y_max)  ; overwrites original y_mid

  title = title+' ['+strtrim(string(n_ss),2)+']'

  plot,[0],[0],$
    /nodata,$
    title=title,$
    xstyle=1,$
    xrange=[min(r),max(r)],$
    xtitle=xtitle,$
    ystyle=1,$
    yrange=y_mid*i_zero+[-i_zero,1]*(y_max-y_mid*i_zero)*zoom,$
    ytitle=ytitle,$
    color=line
; original:    yrange=y_mid*i_zero+[-i_zero,1]*zoom,$

  ;; Equilibrium

  oplot,r,y_eq,color=color_vec[2],psym=8,symsize=0.5

  ;; Equilibrium plus fluctuations

  oplot,r,y0_ave,color=line

  ;; Plot singular surface(s)

  get_singsurf_vec,n_ss,r_surf  
  n_surf = n_elements(r_surf)

  for i=0,n_surf-1 do begin
     oplot,r_surf[i]*[1,1],[-100,100],color=color_vec[3]
  endfor

  if (i_ss_mn_den_plot eq 1) then begin
     get_singsurf_den,r_surf
  endif

  plot_bnd
  plot_finish

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin

     openw,1,pname+'.idlout'

     for i=0,n_r-1 do begin
        printf,1, r[i],y_eq[i],y0_ave[i]
     endfor

     printf,1,'Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
     
  endif
  ;;---------------------------------------------------

end 
