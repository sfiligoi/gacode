pro dist_plot

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_DIST

  t_str = ' t='+strtrim(string(t[t_c]),2)
  k_str = ' [k='+strtrim(string(fix(ik_c)),2)+']'

  title = t_str+k_str

  xtitle = '!4h!3/!4p!3'

  if (i_spec lt n_ion) then begin
     tag_s = 'Ion '+strtrim(string(i_spec+1),2)
     ptag_s = 'ion'+strtrim(string(i_spec+1),2)
  endif else begin
     tag_s = 'Electron'
     ptag_s = 'elec'
  endelse

  ytitle = tag_s+' gyrocenter distribution'

  pname = 'distfn-'+ptag_s

  plot_def_new,pname

  if (ik_c lt n_pass) then begin

     y = fltarr(2,n_theta_p)

     y[*,*] = hp[*,*,ik_c,i_spec,t_c] 

     if (dist_axis eq 0) then begin
        ymax = 1.5*max(abs(y))
        ymin = -ymax
     endif else begin
        ymax = 1.1*max(abs(y))
        ymin = 0.0
     endelse

     plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[-1,1],$
       xtitle=xtitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax],$
       ytitle=ytitle,$
       color=line
     
     oplot,theta_t_p[*,ik_c]/!pi,y[0,*],$
       color=color_vec[0],linestyle=line_vec[0] 
     oplot,theta_t_p[*,ik_c]/!pi,y[0,*],$
       color=pen,psym=8,symsize=symsize 

     oplot,theta_t_p[*,ik_c]/!pi,y[1,*],$
       color=color_vec[1],linestyle=line_vec[1] 
     oplot,theta_t_p[*,ik_c]/!pi,y[1,*],$
       color=pen,psym=8,symsize=symsize 

  endif else begin

     k = ik_c-n_pass

     y = fltarr(2,n_theta_t)

     y[*,*] = ht[*,*,k,i_spec,t_c]

     if (dist_axis eq 0) then begin
        ymax = 1.5*max(abs(y))
        ymin = -ymax
     endif else begin
        ymax = max(y)
        ymin = min(y)
        dy = ymax-ymin
        ymax = ymax+0.2*dy
        ymin = ymin-0.2*dy
     endelse

     plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xminor=0,$
       xrange=[-1,1],$
       xtitle=xtitle,$
       ystyle=1,$
       yminor=0,$
       yrange=[ymin,ymax],$
       ytitle=ytitle,$
       color=line
     
     oplot,theta_t_t[*,k]/!pi,y[0,*],$
       color=color_vec[0],linestyle=line_vec[0] 
     oplot,theta_t_t[*,k]/!pi,y[0,*],$
       color=pen,psym=8,symsize=symsize 

     oplot,theta_t_t[*,k]/!pi,y[1,*],$
       color=color_vec[1],linestyle=line_vec[1] 
     oplot,theta_t_t[*,k]/!pi,y[1,*],$
       color=pen,psym=8,symsize=symsize 

  endelse

  plot_finish

end
