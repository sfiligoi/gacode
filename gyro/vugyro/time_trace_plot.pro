pro time_trace_plot

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_TIME_TRACE
  
  time_phi_labels

  case (index) of 

     1: begin

        pname = 'phi_mag_1'

        plot_def_new,pname
        
        j_str = '!4h/p!3 = '+strtrim(string(thetad_plot[j_c]),2)
        r_str = ' r = '+strtrim(string(r[i_c]),2)
        n_str = ' n = '+strtrim(string(n_tor[in_c]),2)

        title = j_str+r_str+n_str

        ymax = max(abs(u[*,j_c,i_c,i_f,in_c,*]))

        yr = u[0,j_c,i_c,i_f,in_c,*]/ymax
        yi = u[1,j_c,i_c,i_f,in_c,*]/ymax

        plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[min(t),max(t)],$
          xtitle=csa_string,$
          ystyle=1,$
          yminor=0,$
          yrange=[-1,1],$
          ytitle=ytitle[active_f-1],$
          color=line

        oplot,t,yr,color=color_vec[0]   
        oplot,t,yi,color=color_vec[1]  

        oplot,[1,1]*t[1],[1,1]*yr[1],$
          color=1,$
          psym=-8,$
          symsize=dotsize  

        oplot,[1,1]*t[1],[1,1]*yi[1],$
          color=line,$
          psym=-8,$
          symsize=dotsize
        
     end   

     2: begin
        
        plot_def_new,'phi_mag_2'
        
        j_str = '!4h!3 = (all)'
        r_str = ' r = '+strtrim(string(r[i_c]),2)
        n_str = ' n = '+strtrim(string(n_tor[in_c]),2)

        title = j_str+r_str+n_str

        y = fltarr(n_time)

        for tt=0,n_time1 do begin
           
           y[tt] = total(u[0,*,i_c,i_f,in_c,tt]^2 + $
                         u[1,*,i_c,i_f,in_c,tt]^2)

        endfor

        y = sqrt(y)
        ymax = max(y)
        y = y/ymax

        ymax=1.0
        ymin=0.0

        print, 'i_er_shear=',i_er_shear
        
        if (i_er_shear > 0) then begin

           print, 'i_c=',i_c

           for tt=10,n_time1 do begin

              yy0=total(u[0,*,i_c,i_f,in_c,tt]^2 + $
                        u[1,*,i_c,i_f,in_c,tt]^2)
              yyp=total(u[0,*,i_c+1,i_f,in_c,tt]^2 + $
                        u[1,*,i_c+1,i_f,in_c,tt]^2)
              
              yym=total(u[0,*,i_c-1,i_f,in_c,tt]^2 + $
                        u[1,*,i_c-1,i_f,in_c,tt]^2)


              yy0=sqrt(yy0)
              yyp=sqrt(yyp)
              yym=sqrt(yym)

              ;;---------------------------------------------------

              yy0=total(u[0,*,i_c,i_f,in_c,tt])
              yyp=total(u[0,*,i_c+1,i_f,in_c,tt])
              yym=total(u[0,*,i_c-1,i_f,in_c,tt])

              rr0=r[i_c]
              rrp=r[i_c+1]
              rrm=r[i_c-1]

              y[tt] = ((yyp-yy0)/(rrp-rr0)-(yy0-yym)/(rr0-rrm)) $
                /((rrp-rrm)/2.)

              if(i_er_shear eq 2) then begin
                 y[tt]=yy0
              endif

           endfor

           for tt=0,10 do begin
              y[tt]=0.0
           endfor
           
           ymax = max(y)
           ymin = min(y)
           
           i_log=1

        endif

        if (i_log eq 0) then begin

           plot_io,[1],[1],$
             /nodata,$
             title=title,$
             xstyle=1,$
             xminor=0,$
             xrange=[min(t),max(t)],$
             xtitle=csa_string,$
             ystyle=1,$
             yminor=0,$
             yrange=[1e-5,1],$
             ytitle=ytitle[active_f-1],$
             color=line

        endif else begin

           plot,[0],[1],$
             /nodata,$
             title=title,$
             xstyle=1,$
             xminor=0,$
             xrange=[min(t),max(t)],$
             xtitle=csa_string,$
             ystyle=1,$
             yminor=0,$
             yrange=[ymin,ymax],$
             ytitle=ytitle[active_f-1],$
             color=line
           
        endelse

        oplot,t,y,color=color_vec[0]   
        
        oplot,[1,1]*t[t_dot],[1,1]*y[t_dot],$
          color=line,$
          psym=-8,$
          symsize=dotsize

        oplot,[t[0],t[n_time1]],[1,1]*y[t_dot],$
          color=line,$
          linestyle=1,$
          symsize=dotsize
        
     end   

     3: begin

        plot_def_new,'phi_mag_3'
        
        j_str = '!4h/p!3 = '+strtrim(string(thetad_plot[j_c]),2)
        r_str = ' r = (all)'
        n_str = ' n = '+strtrim(string(n_tor[in_c]),2)

        title = j_str+r_str+n_str

        y = fltarr(n_time)

        for tt=0,n_time1 do begin

           y[tt] = total(u[0,j_c,*,i_f,in_c,tt]^2 + $
                         u[1,j_c,*,i_f,in_c,tt]^2)
        endfor

        y = sqrt(y)
        ymax = max(y)
        y = y/ymax

        if (i_log eq 0) then begin

           plot_io,[1],[1],$
             /nodata,$
             title=title,$
             xstyle=1,$
             xminor=0,$
             xrange=[min(t),max(t)],$
             xtitle=csa_string,$
             ystyle=1,$
             yminor=0,$
             yrange=[1e-5,1],$
             ytitle=ytitle[active_f-1],$
             color=line

        endif else begin

           plot,[0],[1],$
             /nodata,$
             title=title,$
             xstyle=1,$
             xminor=0,$
             xrange=[min(t),max(t)],$
             xtitle=csa_string,$
             ystyle=1,$
             yminor=0,$
             yrange=[0,1],$
             ytitle=ytitle[active_f-1],$
             color=line

        endelse

        oplot,t,y,color=color_vec[0]

        oplot,[1,1]*t[t_dot],[1,1]*y[t_dot],$
          color=line,$
          psym=-8,$
          symsize=dotsize

        oplot,[t[0],t[n_time1]],[1,1]*y[t_dot],$
          color=line,$
          linestyle=1,$
          symsize=dotsize

     end

     4: begin

        plot_def_new,'phi_mag_4'
        
        j_str = '!4h!3 = (all)'
        r_str = ' r = (all)'
        n_str = ' n = '+strtrim(string(n_tor[in_c]),2)

        title = j_str+r_str+n_str

        y = fltarr(n_time)

        for tt = 0,n_time1 do begin
           y[tt] = $
             total(U[0,*,*,i_f,in_c,tt]^2+ $
                   U[1,*,*,i_f,in_c,tt]^2)/(n_theta_plot*n_r)
        endfor
        
        y = sqrt(y)
        ymax = max(y)
        ymin = min(y[1:n_time1])

        if (i_log eq 0) then begin

           plot_io,[1],[1],$
             /nodata,$
             title=title,$
             xstyle=1,$
             xminor=0,$
             xrange=[min(t),max(t)],$
             xtitle=csa_string,$
             ystyle=1,$
             yminor=0,$
             yrange=[ymin,ymax],$
             ytitle=ytitle[active_f-1],$
             color=line

        endif else begin

           plot,[0],[1],$
             /nodata,$
             title=title,$
             xstyle=1,$
             xminor=0,$
             xrange=[min(t),max(t)],$
             xtitle=csa_string,$
             ystyle=1,$
             yminor=0,$
             yrange=[0,ymax],$
             ytitle=ytitle[active_f-1],$
             color=line
           
        endelse

        oplot,t,y,color=color_vec[0]   
        
        oplot,[1,1]*t[t_dot],[1,1]*y[t_dot],$
          color=line,$
          psym=-8,$
          symsize=dotsize

        oplot,[t[0],t[n_time1]],[1,1]*y[t_dot],$
          color=line,$
          linestyle=1,$
          symsize=dotsize


     end   

  endcase

  if (index gt 1) then begin
     yr = y[*]
     yi = 0.0*y
  endif

  ;;---------------------------------------------------
  ;; DATA EXPORT
  ;;
  if (plot_export eq 1) then begin
     openw,1,pname+'.idlout'
     for i=0,n_time-1 do begin
        printf,1,t[i],yr[i],yi[i]
     endfor
     close,1
     print,'Exported data to ',pname+'.idlout'
     plot_export = 0
  endif
  ;;---------------------------------------------------

  plot_finish

end
