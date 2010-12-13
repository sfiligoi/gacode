pro midplane_each_power_event, midplane_each_power 

  common GLOBAL
  common PLOT_VARIABLES
  common PRIVATE_MIDPLANE_EACH_POWER
  common MIDPLANE_DATA

  widget_control, midplane_each_power.id, $
    get_uvalue=uvalue

  wset, widget

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 
     
     1: goto, plot_it
     
     6: begin
        bar_plot_n2_flag = 1-bar_plot_n2_flag
        goto, plot_it
     end

     61: begin
        div_dn_flag = 1-div_dn_flag
        goto, plot_it
     end

     7: begin
        plot_mode = 2
        goto, plot_it
     end

     8: begin
        goto, plot_it
     end

     9: widget_control, midplane_each_power.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

   znn = fltarr(n_n)
   xxx = fltarr(n_r)

  znn[*] = 0.0
  xxx[*] = 0.0

  i_pwr_n2=2
  ;; i_pwr = 2 is electron density or ion_2 density
  for i_n=0,n_n-1 do begin
     xxx[*] = 0.0
     for i_time=it1,it2 do begin
        xxx[*] = xxx[*] +abs(fft(pwr[i_pwr_n2,*,i_n,i_time]))^2
     endfor
     xxx[*] = xxx[*]/(it2-it1+1)

     for i=0,n_r-1 do begin
      znn[i_n]=znn[i_n]+xxx[i]
     endfor
  ;;ERROR      znn[i_n]=znn[i_n]/n_r  REW: divide by n_r ?
  endfor

  zn_h = 0.0
  zn_l = znn[0]
  ky_zn_h = 0.0
  ky_zn_l = 0.0


  for i_n=1,n_n-1 do begin
   if (kt_rho[i_n] le  1.0) then begin
    zn_l=zn_l+2.*znn[i_n]
    ky_zn_l=ky_zn_l+2.*znn[i_n]*kt_rho[i_n]
   endif
   if (kt_rho[i_n] gt  1.0) then begin
    zn_h=zn_h+2.*znn[i_n]
    ky_zn_h=ky_zn_h+2.*znn[i_n]*kt_rho[i_n]
   endif
  endfor

   print, '-----------------------------------------------------------------------------'
   print, 'n_r=',n_r,' ','(2*n_n+1)=',2*n_n+1
   print, '-----------------------------------------------------------------------------'
   print, 'rho_star=',rho_s
   print, '-----------------------------------------------------------------------------'

    print, 'kt_rho           znn'
    print,  kt_rho[0], '  ', znn[0]
   for i_n=1,n_n-1 do begin
    print,  kt_rho[i_n], '  ', 2.*znn[i_n]
   endfor
    print, 'zn_l=',zn_l
    print, 'zn_h=',zn_h
   print, '-----------------------------------------------------------------------------'

  yn = fltarr(n_n)
  
  smf_tag,title,pname,ytitle
  title = title+'1D Power Spectrum:  [n_tilda/n0]**2'
  pname = 'n2_n-'+pname

  title = '1D Power Spectrum: [n_tilda/n0]**2/rho_star**2'
  pname = 'n2_n'


     yn[0] = znn[0]/rho_s^2
     for i_n=1,n_n-1 do begin
        yn[i_n] = 2.*znn[i_n]/rho_s^2
     endfor

  if(div_dn_flag eq 1 ) then begin
    for i_n=0,n_n-1 do begin
        yn[i_n] = yn[i_n]/n_dn
     endfor
  endif



  plot_def_new,pname

  x = kt_rho[*]
  xtitle = '!3'+kt_rho_string

     ytitle='!3absolute contribution per mode'

  xmin = min(x)
  xmax = max(x)+0.5*(x[1]-x[0])
  ymin = min(yn)
  ymax = max(yn)*1.3

  plot,[0],[0],$
    /nodata, $
    xtitle=xtitle+' with '+t_string, $
    xminor=1, $
    xstyle=1,xrange=[xmin,xmax],$
    ytitle=ytitle, $
    ystyle=1,yrange=[ymin,ymax],$
    title=title, $
    color=line

  oplot,[xmin,xmax],[0,0],linestyle=1

   if (bar_plot_n2_flag eq 0) then begin
     bar_oplot,x,yn,color_vec[0]
   endif else begin
     oplot,x,yn,psym=8,color=line
     oplot,x,yn,color=color_vec[0]
   endelse

  ;;----------------------------------------------------
  ;; Labels for separate ITG/TEM -ETG_1 -ETG_2  [n_tilda/n0]**2

  n2_itg_tem   = 0.0
  n2_etg_1 = 0.0
  n2_etg_2 = 0.0

  for i_n=0,n_n-1 do begin
     if (x[i_n] gt 2.0) then n2_etg_2 = n2_etg_2+yn[i_n]
     if (x[i_n] gt 1.0) then n2_etg_1 = n2_etg_1+yn[i_n]
     if (x[i_n] le 1.0) then n2_itg_tem = n2_itg_tem+yn[i_n]
  endfor

  print, 'n2_itg_tem=',n2_itg_tem
  print, 'n2_etg_1=',n2_etg_1
  print, 'n2_etg_2=',n2_etg_2

  x0 = xmin + 0.1*(xmax-xmin)

  str = strmid(strtrim(string(n2_itg_tem),2),0,15)
  y0 = ymax - 0.06*(ymax-ymin)
  xyouts,x0,y0,xtitle+' < 1 (ITG/TEM) ['+str+']'

  str = strmid(strtrim(string(n2_etg_1),2),0,15)
  y0 = ymax - 0.11*(ymax-ymin)
  xyouts,x0,y0,xtitle+' > 1 (ETG_1) ['+str+']'

  str = strmid(strtrim(string(n2_etg_2),2),0,15)
  y0 = ymax - 0.16*(ymax-ymin)
  xyouts,x0,y0,xtitle+' > 2 (ETG_2) ['+str+']'
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

  return

end
