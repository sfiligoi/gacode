pro field_r_plot

  common GLOBAL
  common PLOT_VARIABLES

  if (i_abs eq 0) then begin
     prest = ''
  endif else begin
     if (i_log eq 0) then begin
        prest = 'abs '
     endif else begin
        prest = 'log '
     endelse
  endelse

  if (i_f eq 0) then begin
     fst = '!4u!3'
     pname = 'phi_r'
  endif else begin
     fst = '!3A'
     pname = 'ap_r'
  endelse
  
  nst = strtrim(string(fix(n_tor[in_c])),2)
  if i_all eq 0 then begin
     tst = strtrim(string(t_c),2)
  endif else begin
     tst = 'all'
  endelse

  xtitle = '!3r/a'
  title  = prest+fst+': n='+nst+' t='+tst
  ytitle = ' '

  pname = pname+nst+'-t'+tst

  plot_def_new,pname

  z = fltarr(2,n_r,n_time)
  z[*,*,*] = u[*,j_c,*,i_f,in_c,*]

  ;** Add code this way to view other fields.
  ;z[*,*,*] = mom_n[*,0,*,0,in_c,*]

  if i_all eq 1 then begin

     ymax = 1.1*max(abs(z)) 
     ymin = -ymax

     if (i_abs eq 1) then begin
        ymin = 0.0
        if (i_log eq 1) then begin
           ymax = alog(max(abs(z)))
           ymin = alog(min(abs(z)))
        endif
     endif

     plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[min(r),max(r)],$
          xtitle=xtitle,$
          ystyle=1,$
          yminor=0,$
          yrange=[ymin,ymax],$
          ytitle=ytitle,$
          color=line

     for i_time=it1,it2 do begin

        if (i_abs eq 0) then begin
           oplot,r,z[0,*,i_time],color=color_vec[0]
           oplot,r,z[1,*,i_time],color=color_vec[1]
        endif else begin
           u_plot = fltarr(n_r)
           u0_plot = fltarr(n_r)
           u0_plot[*] = u[0,*,i_time]
           u1_plot = fltarr(n_r)
           u1_plot[*] = u[1,*,i_time]
           u_plot = sqrt(u0_plot^2+u1_plot^2)

           if (i_log eq 1) then begin
              u_plot = alog(u_plot)
           endif
           
           oplot,r,u_plot,color=color_vec[0]              
        endelse
        
     endfor

  endif else begin  

     u_plot = fltarr(n_r)
     u0_plot = fltarr(n_r)
     u0_plot[*] = z[0,*,t_c]
     u1_plot = fltarr(n_r)
     u1_plot[*] = z[1,*,t_c]
     u_plot=sqrt(u0_plot^2+u1_plot^2)

     if (i_abs eq 0) then begin
        ymax = 1.1*max(abs(z[*,*,t_c])) 
        ymin=-ymax
     endif

     if (i_abs eq 1) then begin
        ymax = 1.1*max(u_plot)
        ymin = 0.0
        if (i_log eq 1) then begin
           u_plot=alog(u_plot)
           ymax = max(u_plot)
           ymin = min(u_plot)
        endif
     endif

     plot,[0],[0],$
          /nodata,$
          title=title,$
          xstyle=1,$
          xminor=0,$
          xrange=[min(r),max(r)],$
          xtitle=xtitle,$
          ystyle=1,$
          yminor=0,$
          yrange=[ymin,ymax],$
          ytitle=ytitle,$
          color=line

     if (i_abs eq 0) then begin
        oplot,r,u0_plot,color=color_vec[0]
        oplot,r,u1_plot,color=color_vec[1]
        if (plot_mode eq 1) then begin
           oplot,r,u0_plot,color=line,psym=8,symsize=dotsize
           oplot,r,u1_plot,color=line,psym=8,symsize=dotsize
        endif
     endif else begin
        oplot,r,u_plot,color=color_vec[0] 
     endelse

     ;;---------------------------------------------------
     ;; DATA EXPORT
     ;;
     if (plot_export eq 1) then begin
        openw,1,pname+'.idlout'
        for i=0,n_r-1 do begin
           printf,1,r[i],z[*,i,t_c]
        endfor
        close,1
        print,'Exported data to ',pname+'.idlout'
        plot_export = 0
     endif
     ;;---------------------------------------------------

  endelse

  plot_finish

  return

end
