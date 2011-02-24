pro gbflux_i_plot

  common GLOBAL
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA
  common ZMOMENT_DATA
  common PRIVATE_GBFLUX_I

  ;;-------------------------------------------------------
  ;; Text labels
  ;;
  gbflux_tag,title,pname,ytitle,units1,units2,units3
  
  pname = 'gbflux_i-'+pname
  
  case (i_units) of

     0: begin
        ;; flux in gyrobohm units
        xnorm = 1.0
        units = units1
     end

     1: begin
        ;; flux in W/m^2, etc.
        xnorm = xunits[9+i_moment]
        units = units2
     end 

     2: begin 
        ;; flux in W, etc.
        xnorm = xunits[9+i_moment]
        units = units3
     end

  endcase

  xtitle = '!3r/a with'+t_string

  plot_def_new,pname
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Pick out correct diffusivity and get the time 
  ;; average.
  ;;
  y_i = fltarr(n_r,n_time)

  if i_f ge 0 then begin
     y_i[*,*] = gbflux_i[i_spec,i_f,i_moment,*,*]
  endif else begin
     y_i[*,*] = 0.0
     for ix=0,n_field-1 do begin
        y_i[*,*] = y_i[*,*]+gbflux_i[i_spec,ix,i_moment,*,*]
     endfor
  endelse

  y   = fltarr(n_time)
  y_r = fltarr(n_r)

  for i=0,n_r-1 do begin
     y[*] = y_i[i,*]
     diff_stat_fast,y,it1,it2,ave_y
     y_r[i] = ave_y     
  endfor
  ;;-------------------------------------------------------
  
  ;;-------------------------------------------------------
  ;; Set the appropriate experimental chi-profile
  ;;
  if (exists_exp_derived) then begin

     y_exp = fltarr(n_rho)
     y_exp[*] = 0.0

     if (i_moment eq 0) then begin

        ;; Particle flux

        y_exp[*] = exp_profile[15,*]/xunits[9]/exp_derived[22,*]

     endif

     if (i_moment eq 1) then begin

        ;; Energy flux

        if (i_spec eq n_kinetic-1) or (electron_method eq 3) then begin
           ;; Qe
           y_exp[*] = exp_profile[11,*]/xunits[10]/exp_derived[22,*]
        endif else begin
           ;; Qi
           y_exp[*] = exp_profile[12,*]/xunits[10]/exp_derived[22,*]
        endelse
        
     endif

     if (i_moment eq 2) then begin

        ;; Momentum flux

        y_exp[*] = exp_profile[10,*]/xunits[11]/exp_derived[22,*]
        
     endif

  endif

  ;; ** At this point, y and y_exp are in gyroBohm units **

  if (i_units eq 2) then begin
     ;; Multiply by V'
     y_r[*]   = y_r[*]*INTERPOL(exp_derived[22,*],r_from_rho,r)
     y_exp[*] = y_exp[*]*exp_derived[22,*]
  endif
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
       yrange=[-i_zero*zoom,zoom]*max(abs(y_r*xnorm)),$
       ytitle=ytitle+units,$
       color=line

  oplot,r,y_r*xnorm,color=color_vec[0]
  if (exists_exp_derived) then oplot,r_from_rho,y_exp*xnorm,color=color_vec[1]

  ;; Ron's addition
  if (i_moment eq 3) then begin
     if (i_units eq 2) then begin
        y_r2 = fltarr(n_r)
        y_r2[*] = 0.0

        ;; y_r2 is the radial integral over simulation  MW/m --->MW
        for i=8,n_r-8 do begin
           y_r2[i]= y_r2[i-1]+(y_r[i]+y_r[i-1])/2.*xunits[2]*(r[n_r-1]-r[0])/n_r
        endfor
        print, 'plot radially integrated exchange flow in MW over simulation'
        oplot,r,y_r2*xnorm,color=color_vec[1]
     endif
  endif


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
     y_exp_fine = INTERPOL(y_exp*xnorm,r_from_rho,r)
     openw,1,pname+'.idlout'
     for i=0,n_r-1 do begin
        printf,1,r[i],y_r[i]*xnorm,y_exp_fine[i]
     endfor
     printf,1,'Average:',t_min,t_max
     printf,1,'Units: '+units
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  return

end
