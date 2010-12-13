pro balloon_plot,i_field

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_BALLOON
  common PROFILE_SIM_DATA

  ytitle  = balloon_tag[i_field]
  tfile   = balloon_tag[i_field]
  y_factor = 1.1*dy_balloon

  title = '!3k!4!dh!nq!d!3s!n!4 = '+$
          strmid(strtrim(string(kt_rho[0]),2),0,5)+$
          '!3  l = '+strtrim(string(balloon_l0),2)

  plot_def_new,tfile

  xpmax = 2*i_p+1 
  xpmin = -xpmax 

  if (i_p eq n_r/2) then xpmax = xpmax-2

  plot,[0],[0],$
       /nodata,$
       title=title,$
       xstyle=1,$
       xtitle='!4h!d*!n/p!3',$
       xrange=[xpmin,xpmax],$
       ystyle=1,$
       yrange=[-1,1]*y_factor,$
       ytitle=ytitle,$
       color=line
  
  if balloon_norm eq 0 then begin

     ;; Original GYRO normalization

     c_r = 1.0

  endif else begin

     ;; Normalize to selected eigenmode

     c_r = max(abs(balloon_plot[*,balloon_next/2,0,balloon_norm-1,t_c]))

  endelse

  x = -(1+balloon_np)+2*balloon_np*findgen(balloon_next)/balloon_next

  oplot,x,balloon_plot[0,*,balloon_l0,i_field,t_c]/c_r,$
        color=color_vec[1],$
        linestyle=0

  oplot,x,balloon_plot[1,*,balloon_l0,i_field,t_c]/c_r,$
        color=color_vec[2],$
        linestyle=0

  oplot,x,balloon_plot[0,*,balloon_l0,i_field,t_c]/c_r,$
        color=line,$
        psym=8,$ 
        symsize=dotsize

  oplot,x,balloon_plot[1,*,balloon_l0,i_field,t_c]/c_r,$
        color=line,$
        psym=8,$ 
        symsize=dotsize

  if (plot_export eq 1) then begin

     openw,1,tfile+'.idlout'

     for j=0,balloon_next-1 do begin
        printf,1,x[j],$
               balloon_plot[0,j,balloon_l0,i_field,t_c]/c_r,$
               balloon_plot[1,j,balloon_l0,i_field,t_c]/c_r
     endfor

     close,1
     print,'Exported data to ',tfile+'.idlout',t_c
     plot_export = 0

  endif

  for p=-1-balloon_np,balloon_np,2 do begin
     oplot,p*[1,1],[-1,1],color=line,linestyle=1
  endfor

  plot_finish

end
