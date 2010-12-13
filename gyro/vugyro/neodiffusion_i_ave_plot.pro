pro neodiffusion_i_ave_plot

  common GLOBAL
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA
  common ZMOMENT_DATA
  common PRIVATE_NEODIFFUSION_I_AVE

  y_i = fltarr(n_r,n_time)
  y_i[*,*] = 0.0

  if (i_turb_neo ne 5) and (i_moment eq 2) then i_moment=0

  case (i_turb_neo) of

     0: begin

        smf_tag,title,pname,ytitle

        title = title+' diffusion'
        pname = 'diff_i-'+pname

        if i_f ge 0 then begin
           y_i[*,*] = diff_i[i_spec,i_f,i_moment,*,*]
        endif else begin
           y_i[*,*] = diff_i[i_spec,0,i_moment,*,*]+$
             diff_i[i_spec,1,i_moment,*,*]
        endelse

     end

     1: begin

        title = 'Turb. Ion+Electron Energy Diffusivity/2'
        pname = 'diff_i_turb_ie_energy'
        ytitle = 'Diffusivity'
        if i_f ge 0 then begin
           y_i[*,*] = 0.0
           for ii_spec = 0,n_kinetic-1 do begin
              y_i[*,*] = y_i[*,*]+diff_i[ii_spec,i_f,1,*,*]
           endfor
        endif else begin
           y_i[*,*] = 0.0
           for ii_spec = 0,n_kinetic-1 do begin
              y_i[*,*] = y_i[*,*]+diff_i[ii_spec,0,1,*,*]+$
                diff_i[ii_spec,1,1,*,*]
           endfor
        endelse
        y_i[*,*] =  0.5*y_i[*,*]

     end

     2: begin

        title = 'Turb+Neo Ion Energy Diffusivity'
        pname = 'diff_i_turb_neo_i_energy'
        ytitle='Diffusivity'
        if i_f ge 0 then begin
           y_i[*,*] = 0.0
           for ii_spec = 0,n_ion-1 do begin
              y_i[*,*] = y_i[*,*]+diff_i[ii_spec,i_f,1,*,*]
           endfor
        endif else begin
           y_i[*,*] = 0.0
           for ii_spec = 0,n_ion-1 do begin
              y_i[*,*] = y_i[*,*]+diff_i[ii_spec,0,1,*,*]+$
                diff_i[ii_spec,1,1,*,*]
           endfor
        endelse
        if (exists_diff_neo_i eq 1) then begin
           y_i[*,*] = y_i[*,*]+diff_neo_i[0,1,*,*]
        endif

     end

     3: begin 
        title = 'Turb+Neo Electron Energy Diffusivity'
        pname = 'diff_i_turb_neo_e_energy'
        ytitle='Diffusivity'
        if i_f ge 0 then begin
           y_i[*,*] = 0.0
           y_i[*,*] = y_i[*,*]+diff_i[n_kinetic-1,i_f,1,*,*]
        endif else begin
           y_i[*,*] = 0.0
           y_i[*,*] = y_i[*,*]+diff_i[n_kinetic-1,0,1,*,*]+$
             diff_i[i_spec,1,1,*,*]
        endelse
        if (exists_diff_neo_i eq 1) then begin
           y_i[*,*] = y_i[*,*]+diff_neo_i[n_kinetic-1,1,*,*]
        endif

     end

     4: begin

        title = 'Total Energy Diffusivity'
        pname = 'diff_i_turb_ie_energy'
        ytitle = 'Diffusivity'
        if i_f ge 0 then begin
           y_i[*,*] = 0.0
           for ii_spec = 0,n_kinetic-1 do begin
              y_i[*,*] = y_i[*,*]+diff_i[ii_spec,i_f,1,*,*]
           endfor
        endif else begin
           y_i[*,*] = 0.0
           for ii_spec = 0,n_kinetic-1 do begin
              y_i[*,*] = y_i[*,*]+diff_i[ii_spec,0,1,*,*]+$
                diff_i[ii_spec,1,1,*,*]
           endfor
        endelse
        if (exists_diff_neo_i eq 1) then begin
           y_i[*,*] = y_i[*,*]+diff_neo_i[0,1,*,*]
        endif

     end

     5: begin

        smf_tag,title,pname,ytitle

        title = title+' Diffusion_neo'
        pname = 'diff_neo_i-'+pname

        if (i_moment eq 2) then title = 'parallel velocity/cs'

        if (exists_diff_neo_i eq 1) then begin
           y_i[*,*] = diff_neo_i[i_spec,i_moment,*,*]
        endif

     end

  endcase

  plot_def_new,pname

  if (i_units eq 0) then begin
     plot_units = 1.0
     units = ' [units of (c!ds!n/a)!4q!3!ds!n!u2!n]'
  endif else begin
     plot_units = xunits[8]
     units = ' [units of (m!u2!n/sec)]'
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

  if (i_div eq 1) then begin

     gradient_factor = fltarr(n_r)

     y  = fltarr(n_r,n_time)
     y0 = fltarr(n_r)
     y1 = fltarr(n_r,n_time)

     if (i_moment eq 0) then begin
        y0[*] = dlnndr_s[i_spec,*]
        y1[*,*] = zmoment[*,i_spec,0,*]
        vector_deriv,r,y1,y,boundary_method
        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+den_s[i_spec,i]*y0[i]
           y[i,*] = y[i,*]/(den_s[i_spec,i]+y1[i])
        endfor
     endif else begin
        y0[*] = dlntdr_s[i_spec,*]
        for i=0,n_r-1 do begin
           y1[i,*] = (zmoment[i,i_spec,1,*]- $
                      1.5*zmoment[i,i_spec,0,*]*tem_s[i_spec,i])/ $
             (1.5*(den_s[i_spec,i]+zmoment[i,i_spec,0,*]))        
        endfor
        vector_deriv,r,y1,y,boundary_method
        for i=0,n_r-1 do begin
           y[i,*] = -y[i,*]+tem_s[i_spec,i]*y0[i]
           y[i,*] = y[i,*]/(tem_s[i_spec,i]+y1[i])
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
  
  if (i_div eq 2) then begin

     div_flow = fltarr(n_r)
     div_flow[*] = 0.0

     y_s = fltarr(n_r)
     y_s[*] = y_r[*]
     
     ;; smoother
     for i=1,n_r-2 do begin
        y_s[i]= y_r[i]/2.0+(y_r[i+1]+ y_r[i-1])/4.
     endfor

     ntdlntdr_s =  fltarr(n_r)
     for i=n_bnd,n_r-1-n_bnd do begin
        ntdlntdr_s[i]= den_s[i_spec,i]*dlntdr_s[i_spec,i]*tem_s[i_spec,i]
     endfor

     ntdlnndr_s =  fltarr(n_r)
     for i=n_bnd,n_r-1-n_bnd do begin
        ntdlnndr_s[i]= den_s[i_spec,i]*dlnndr_s[i_spec,i]
     endfor

     ntdlntndr_s =  fltarr(n_r)
     if (i_moment eq 0) then begin
        ntdlntndr_s[*] = ntdlnndr_s[*]
     endif
     if (i_moment eq 1) then begin
        ntdlntndr_s[*] = ntdlntdr_s[*]
     endif
     
     for i=n_bnd,n_r-1-n_bnd do begin
        div_flow[i] = r[i+1]*y_s[i+1]*ntdlntndr_s[i+1]
        div_flow[i] = div_flow[i]-r[i-1]*y_s[i-1]*ntdlntndr_s[i-1]
        div_flow[i] = div_flow[i]/r[i]/(r[i+1]-r[i-1]) 
     endfor

  endif

  xtitle = '!3r/a with'+t_string
  title = title+' ['+strtrim(string(n_ss),2)+']'

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
  
  if i_div eq 2 then begin
     oplot,r,div_flow,color=color_vec[1]
  endif

  case (i_turb_neo) of 

     0: begin
        y_exp[*]=0.0
        if(i_moment eq 1) then begin
           if (i_spec eq n_kinetic-1) then begin
              y_exp[*] = chi_e_exp[*]*plot_units
           endif else begin
              y_exp[*] = chi_i_exp[*]*plot_units
           endelse
        endif
        if(i_moment eq 0) then begin
           if (i_spec eq n_kinetic-1) then begin
              y_exp[*] = diff_ne_exp[*]*plot_units
           endif
        endif

     end

     1: y_exp[*] = 0.5*(chi_e_exp[*]+chi_i_exp[*])*plot_units

     2: y_exp[*] = chi_i_exp[*]*plot_units

     3: y_exp[*] = chi_e_exp[*]*plot_units

     4: y_exp[*] = chi_e_exp[*]*plot_units

     5: y_exp[*] = 0.0

  endcase

  oplot,r,y_exp,color=color_vec[1]

  ;; Plot singular surface(s)

  get_singsurf_vec,n_ss,r_surf  
  n_surf = n_elements(r_surf)

  for i=0,n_surf-1 do begin
     oplot,r_surf[i]*[1,1],[-100,100],color=color_vec[3]
  endfor

  plot_bnd

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
