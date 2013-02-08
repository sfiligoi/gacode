; Makes
; gbflux-elec-density-phi.ps
; gbflux-elec-energy-phi.ps
; gbflux-ion1-density-phi.ps
; gbflux-ion1-energy-phi.ps
; gbflux-ion2-density-phi.ps
; gbflux-ion2-energy-phi.ps
;  ... for all ions

pro gbflux_auto

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_GBFLUX  

; Initialize variables that control plots
i_f=0       ; up to n_field
i_zero=0    ; 0 or 1
zoom=1.
i_units=0   ; 0, 1, 2
i_tp=0   ; 0:all, 1:trapped, 2:passing
n_ss=0


;  Do particle fluxes first, then energy fluxes:
; i_moment=0  ; up to p_moment-1; density, energy, momentum, exchange
; i_spec=0    ; up to n_kinetic-1, electrons are last

for i_moment=0,1 do begin
  for i_spec=0,n_kinetic-1 do begin 
; plot_mode: anything but 1 produces a .ps file
     plot_mode=2

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  y = fltarr(n_time)
  
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

  case (i_tp) of

     0: begin

        pname = 'gbflux-'+pname
        title = title
        if i_f ge 0 then begin
           y[*] = gbflux[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = 0.0
           for ix=0,n_field-1 do begin
              y[*] = y[*]+gbflux[i_spec,ix,i_moment,*]
           endfor
        endelse

     end

     1: begin

        pname = 'gbflux_t-'+pname
        title = title+' (trapped)'
        if i_f ge 0 then begin
           y[*] = gbflux_trapped[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = 0.0
           for ix=0,n_field-1 do begin 
              y[*] = y[*]+gbflux_trapped[i_spec,ix,i_moment,*]
           endfor
        endelse

     end

     2: begin

        pname = 'gbflux_p-'+pname
        title = title+' (passing)'
        if i_f ge 0 then begin
           y[*] = gbflux[i_spec,i_f,i_moment,*]-$
                  gbflux_trapped[i_spec,i_f,i_moment,*]
        endif else begin
           y[*] = 0.0
           for ix=0,n_field-1 do begin
              y[*] = y[*] + gbflux[i_spec,ix,i_moment,*]-$
                     gbflux_trapped[i_spec,ix,i_moment,*] 
           endfor
        endelse

     end

  endcase

  plot_def_new,pname

  ;; Set units here

  y[*] = y[*]*xnorm
  if (i_units eq 2) then begin
     ;; Multiply by V'- CH fix 12.20.2012
;     y[*] = y[*]*INTERPOL(exp_derived[22,*],r_from_rho,r[n_r/2])
     y[*] = y[*]*INTERPOL(exp_derived[23,*],r_from_rho,r[n_r/2])
  endif

  y_axis_magic,y,ymin,ymax,d_y

  unity = fltarr(it2-it1)
  unity[*] = 1.0

  if (it2-it1) gt 20 then begin

     ;; Enough points for Mikkelsen's method

     res = meanstddev(pname,0,y,t,t_min,t_max)
     ave_y   = res.avg
     err_y   = res.sig2
     trend_y = res.trend
     hump_y  = res.hump

  endif else begin

     ;; Only compute simple average

     diff_stat_fast,y,it1,it2,ave_y
     err_y   = -1
     trend_y = -1
     hump_y  = -1

  endelse

  ave_str   = strmid(strtrim(string(ave_y),2),0,5)
  err_str   = strmid(strtrim(string(err_y),2),0,5)
  trend_str = strmid(strtrim(string(trend_y),2),0,5)
  hump_str  = strmid(strtrim(string(hump_y),2),0,5)

  text = ' !3['+ave_str+'!9+!3'+err_str+']'

  plot,[0],[0],$
       /nodata,$
       title=title+text,$
       xstyle=1,$
       xminor=0,$
       xrange=[min(t),max(t)],$
       xtitle=csa_string,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax]/zoom,$
       ytitle=ytitle+units, $
       color=line
  
  ;; flux trace

  oplot,t,y,color=color_vec[0]

  ;; Average line

  oplot,t[it1:it2],unity*ave_y,color=color_vec[1]   

  ;; Curve dots

  oplot,t[it1]*[1,1],y[it1]*[1,1],psym=8,color=line   
  oplot,t[it2]*[1,1],y[it2]*[1,1],psym=8,color=line   

  ;; RMS deviation bars

  oplot,t[it1:it2],unity*(ave_y-err_y),color=color_vec[2]   
  oplot,t[it1:it2],unity*(ave_y+err_y),color=color_vec[2]   

  dy = ymax-ymin

  xyouts,0.7*max(t),(ymin+0.95*dy)/zoom,$
         '!3trend='+trend_str

  xyouts,0.7*max(t),(ymin+0.9*dy)/zoom,$
         '!3hump='+hump_str

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     printf,1,ave_y,err_y
     printf,1,' Average:',t_min,t_max
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
