pro midplane_power_plot

  common GLOBAL
  common MIDPLANE_DATA
  common PLOT_VARIABLES
  

  if r_flag eq 0 then begin

     ;;=============================
     ;; k_theta rho_s power spectrum
     ;;=============================
     
     y  = fltarr(n_n)
     yt = fltarr(n_time)

     title = tag[i_pwr]+' Power Spectrum'
     pname = 'pwr_n_allr_'+tag[i_pwr]

     for i_n=0,n_n-1 do begin
        yt[*] = 0.0
        for i=0,n_r-1 do begin
           yt[*] = yt[*]+abs(pwr[i_pwr,i,i_n,*])^2
        endfor
        diff_stat_fast,yt,it1,it2,ave_y
        y[i_n] = ave_y/n_r
     endfor

     plot_def_new,pname
     ytitle = '!310!u4!n <f!u2!n>'

     ;; SCALE 
     y = 1e4*y

     x = kt_rho[*]
     xtitle = '!3'+kt_rho_string+' with'+t_string

     if i_loglog eq 0 then begin

        plot,[0],[0],$
          /nodata, $
          xtitle=xtitle, $
          xminor=1, $
          xstyle=1,xrange=[0,max(x)],$
          ytitle=ytitle, $
          ystyle=1,yrange=[0,max(y)],$
          title=title, $
          color=line

        oplot,x,y,psym=8,color=line
        oplot,x,y,color=color_vec[0]

     endif else begin

        plot_io,[1],[1],$
          /nodata, $
          xtitle=xtitle, $
          xminor=1, $
          xstyle=1,xrange=[0,max(x)],$
          ytitle=ytitle, $
          ystyle=1,yrange=[min(y),max(y)],$
          title=title, $
          color=line

        oplot,x,y,psym=8,color=line
        oplot,x,y,color=color_vec[0]

     endelse

     ;;---------------------------------------------------
     ;; DATA EXPORT
     ;;
     if (plot_export eq 1) then begin
        openw,1,pname+'.idlout'
        for i_n=0,n_n-1 do begin
           printf,1,x[i_n],y[i_n]
        endfor
        printf,1,'Average:',t_min,t_max
        close,1
        print,'Exported data to ',pname+'.idlout'
        plot_export = 0
     endif
     ;;---------------------------------------------------

  endif else begin

     ;;=============================
     ;; radial power spectrum
     ;;=============================
     
     n_r_plot = n_r

     y = fltarr(n_r)
     yrt = fltarr(n_r,n_time)     
     yt = fltarr(n_time)     

     title = tag[i_pwr]+' Power Spectrum'
     pname = 'pwr_r_alln_'+tag[i_pwr]

     for i_time=0,n_time-1 do begin
        yrt[*,i_time] = abs(fft(pwr[i_pwr,*,0,i_time]))^2
        for i_n=1,n_n-1 do begin
           yrt[*,i_time] = yrt[*,i_time]+$
             2*abs(fft(pwr[i_pwr,*,i_n,i_time]))^2
        endfor
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

     if (i_loglog eq 0) then begin

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

     endif else begin

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

        endelse

     endelse

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

  endelse

  plot_finish

  return

end
