; Make both .ps and .idlout
; elec-density-phi.ps
; elec-energy-phi.ps
; ion1-density-phi.ps
; ion1-energy-phi.ps
; ... for all ions

pro gbflux_n_auto

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_GBFLUX_N

; Initialize variables that control plots
i_f=0       ; up to n_field
i_units=0   ; 0, 1, 2

frac_flag=0     ; 1:fraction of total flux,   else:not the fraction
bar_plot_flag=1 ; 0:bar,  else:line


;  Do particle fluxes first, then energy fluxes:
i_moment=0  ; up to p_moment-1; density, energy, momentum, exchange
i_spec=0    ; up to n_kinetic-1, electrons are last

for i_moment=0,1 do begin

  for i_spec=0,n_kinetic-1 do begin 
;   Stack IDL commands with &
; plot_mode: anything but 1 produces a .ps file
; Write .idlout file with plot_export=1; can also make .ps file
    plot_mode=2 & plot_export=1


  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  yn = fltarr(n_n)
  y  = fltarr(n_time)
  
  gbflux_tag,title,pname,ytitle,units1,units2,units3

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

  if i_f ge 0 then begin

         for i_n=0,n_n-1 do begin
           y[*] = gbflux_n[i_spec,i_f,i_moment,i_n,*]
           diff_stat_fast,y,it1,it2,y_ave
           yn[i_n] = y_ave
        endfor 

  endif else begin

       for i_n=0,n_n-1 do begin
           y[*] = 0.0
           for ix=0,n_field-1 do begin
              y[*] = y[*]+gbflux_n[i_spec,ix,i_moment,i_n,*]
           endfor
           diff_stat_fast,y,it1,it2,y_ave
           yn[i_n] = y_ave
        endfor

  endelse

  ;; Set units here

  yn[*] = yn[*]*xnorm
  if (i_units eq 2) then begin
     ;; Multiply by V'- CH fix 12.20.12
;     yn[*] = yn[*]*INTERPOL(exp_derived[22,*],r_from_rho,r[n_r/2])
     yn[*] = yn[*]*INTERPOL(exp_derived[23,*],r_from_rho,r[n_r/2])
  endif

  plot_def_new,pname

  x = kt_rho[*]
  xtitle = '!3'+kt_rho_string

  if frac_flag eq 1 then begin
     ytitle='!3fractional contribution per mode'
     yn[*] = yn[*]/total(yn)
  endif

  xmin = min(x)
  xmax = max(x)+0.5*(x[1]-x[0])
  ymin = min(yn)
  ymax = max(yn)*1.3
  if min(yn) gt 0.0 then ymin = 0.0
  if min(yn) lt 0.0 then ymin = min(yn)*1.3

  plot,[0],[0],$
    /nodata, $
    xtitle=xtitle+' with '+t_string, $
    xminor=1, $
    xstyle=1,xrange=[xmin,xmax],$
    ytitle=ytitle+units, $
    ystyle=1,yrange=[ymin,ymax],$
    title=title, $
    color=line

  oplot,[xmin,xmax],[0,0],linestyle=1

  if (bar_plot_flag eq 0) then begin
     bar_oplot,x,yn,color_vec[0]
  endif else begin
     oplot,x,yn,psym=8,color=line
     oplot,x,yn,color=color_vec[0]
  endelse

  ;;----------------------------------------------------
  ;; Labels for separate ITG-ETG diffusivity

  d_itg   = 0.0
  d_etg_1 = 0.0
  d_etg_2 = 0.0
  for i_n=0,n_n-1 do begin
     if (x[i_n] gt 2.0) then d_etg_2 = d_etg_2+yn[i_n]
     if (x[i_n] gt 1.0) then d_etg_1 = d_etg_1+yn[i_n]
     if (x[i_n] le 1.0) then d_itg = d_itg+yn[i_n]
  endfor

  x0 = xmin + 0.1*(xmax-xmin)

  str = strmid(strtrim(string(d_itg),2),0,5)
  y0 = ymax - 0.06*(ymax-ymin)
  xyouts,x0,y0,xtitle+' < 1 (ITG) ['+str+']'

  str = strmid(strtrim(string(d_etg_1),2),0,5)
  y0 = ymax - 0.11*(ymax-ymin)
  xyouts,x0,y0,xtitle+' > 1 (ETG) ['+str+']'

  str = strmid(strtrim(string(d_etg_2),2),0,5)
  y0 = ymax - 0.16*(ymax-ymin)
  xyouts,x0,y0,xtitle+' > 2 (ETG) ['+str+']'
  ;;----------------------------------------------------

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i_n=0,n_n-1 do begin
        printf,1,x[i_n],yn[i_n]
     endfor
     printf,1,'Average:',t_min,t_max
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

  endfor ; loop over species
endfor ; loop over moments

  return

end
