pro midplane_n_power_plot

  common GLOBAL
  common MIDPLANE_DATA
  common PLOT_VARIABLES
  common PROFILE_SIM_DATA
  
  ;;=============================
  ;; radial power spectrum
  ;;=============================

  n_r_plot = n_r
  
  y   = fltarr(n_r)
  yrt = fltarr(n_r,n_time)     
  yt  = fltarr(n_time)     

  n_str = strtrim(string(n_tor[in_c]),2)

  title = tag[i_pwr]+' Power Spectrum (n='+n_str+')'
  pname = 'pwr_r_n'+n_str+'_'+tag[i_pwr]

  for i_time=0,n_time-1 do begin
     yrt[*,i_time] = abs(fft(pwr[i_pwr,*,in_c,i_time]))^2
  endfor

  for i=0,n_r-1 do begin
     yt[*] = yrt[i,*]
     diff_stat_fast,yt,it1,it2,ave_y
     y[i] = ave_y
  endfor

  plot_def_new,pname
  ytitle = '!310!u4!n <f!u2!n>'

  ;; SCALE 
  y = 1e4*shift(y,n_r/2-1)
  x = kr_rho

  xtitle = kr_rho_string+' with'+t_string

  ;; Rational surface harmonics
  kr_rat = fltarr(n_n)
  kr_rat[*] = kt_rho[*]*2*!Pi*shat_s[n_r/2]

  case (i_loglog) of

     0: begin

        ;; LINEAR

        plot,[0],[0],$
             /nodata, $
             xtitle=xtitle, $
             xminor=1, $
             xstyle=1,xrange=[-max(abs(x)),max(abs(x))],$
             ytitle=ytitle, $
             ystyle=1,yrange=[0,max(y)],$
             title=title, $
             color=line

        oplot,x,y,psym=8,color=line
        oplot,x,y,color=color_vec[0]

     end

     1: begin

        if (i_abs eq 0) then begin

           ;; LOG

           plot_io,[1],[1],$
                   /nodata, $
                   xtitle=xtitle, $
                   xminor=1, $
                   xstyle=1,xrange=[-max(abs(x)),max(abs(x))],$
                   ytitle=ytitle, $
                   ystyle=1,yrange=[min(y),max(y)],$
                   title=title, $
                   color=line

           oplot,x,y,psym=8,color=line
           oplot,x,y,color=color_vec[0]

        endif else begin

           ;; LOG-LOG

           xl = x[n_r/2:n_r-1]
           yl = y[n_r/2:n_r-1]
           x = fltarr(n_r/2)        
           y = fltarr(n_r/2)        
           x = xl
           y = yl

           n_r_plot = n_r/2

           plot,[1],[1],$
                /xlog,/ylog,$
                /nodata, $
                xtitle=xtitle, $
                xminor=1, $
                xstyle=1,xrange=[min(x),max(x)],$
                ytitle=ytitle, $
                ystyle=1,yrange=[min(y),max(y)],$
                title=title, $
                color=line

           oplot,x,y,psym=8,color=line
           oplot,x,y,color=color_vec[0]
           for i=1,n_n-1 do begin           
              oplot,kr_rat[i]*[1,1],[1e-10,1],color=color_vec[1]
           endfor

        endelse

     end

  endcase

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_r_plot-1 do begin
        printf,1,x[i],y[i]
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
