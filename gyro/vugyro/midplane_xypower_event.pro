pro midplane_xypower_event, midplane_xypower

  common GLOBAL
  common MIDPLANE_DATA

  widget_control, midplane_xypower.id, $
  get_uvalue=uvalue

  wset, midplane_xypower_wid

  ;;-------------------------------------------------------
  ;; MENU
  ;;-------------------------------------------------------

  case (uvalue) of 

     1: begin
        goto, plot_it
     end   

     2: begin
        counter_up,i_pwr,n_pwr-1,1
        goto, plot_it
     end   

     3: begin
        counter_dn,i_pwr,0,1
        goto, plot_it
     end   

     4: begin
        quad_lin_flag = 1-quad_lin_flag
        goto, plot_it
     end

     5: begin
        log_flag = 1-log_flag
        goto, plot_it
     end

     6: begin
        plot_mode = 2
        goto, plot_it
     end

     8: widget_control, midplane_xypower.top, /destroy

  endcase

  return

  plot_it:

  ;;-------------------------------------------------------
  ;; PLOTTING
  ;;-------------------------------------------------------

  z = fltarr(n_r,n_n)     
  x = fltarr(n_r)
  y = fltarr(n_n)

  if (log_flag eq 0) then begin
     logstr = ' linear'
  endif else begin
     logstr = ' log' 
  endelse
  if (quad_lin_flag eq 0) then begin
     quadlinstr = ' quadratic'
  endif else begin
     quadlinstr = ' linear'
  endelse

  znn = fltarr(n_n)
  xxx = fltarr(n_r)

  znn[*] = 0.0
  xxx[*] = 0.0
  
  for i_n=0,n_n-1 do begin
     xxx[*] = 0.0
     for i_time=it1,it2 do begin
        xxx[*] = xxx[*] +abs(fft(pwr[i_pwr,*,i_n,i_time]))^2
     endfor
     xxx[*] = xxx[*]/(it2-it1+1)

     for i=0,n_r-1 do begin
        znn[i_n]=znn[i_n]+xxx[i]
     endfor
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
  print, 'quadratic spec power low  k*rho_si <= 1 =', zn_l, '%=',zn_l/(zn_l+zn_h)
  print, 'quadratic spec power high k*rho_si >  1 =', zn_h, '%=',zn_h/(zn_l+zn_h)
  print, '-----------------------------------------------------------------------------'
  print, 'quadratic spec ky*power low  k*rho_si <= 1 =', ky_zn_l, '%=',ky_zn_l/(ky_zn_l+ky_zn_h)
  print, 'quadratic spec ky*power high k*rho_si >  1 =', ky_zn_h, '%=',ky_zn_h/(ky_zn_l+ky_zn_h)
  print, '-----------------------------------------------------------------------------'
  print, '-----------------------------------------------------------------------------'
  print, 'n_r=',n_r,' ','(2*n_n+1)=',2*n_n+1,' ','i_pwr=',i_pwr
  print, '-----------------------------------------------------------------------------'
  
  print, 'kt_rho           znn'
  print,  kt_rho[0], '  ', znn[0]
  for i_n=1,n_n-1 do begin
     print,  kt_rho[i_n], '  ', 2.*znn[i_n]   
  endfor
  print, 'zn_l=',zn_l
  print, 'zn_h=',zn_h

  

  title = tag[i_pwr]+quadlinstr+logstr+' 2D Power Spectrum'
  pname = 'pwr_xy_'+tag[i_pwr]

  z[*,*] = 0.0
  for i_n=0,n_n-1 do begin 
     for i_time=it1,it2 do begin       
        z[*,i_n] = z[*,i_n]+abs(fft(pwr[i_pwr,*,i_n,i_time]))^2
     endfor
     z[*,i_n] = shift(z[*,i_n],n_r/2-1)/(it2-it1+1)

     if(quad_lin_flag eq 1) then z[*,i_n] = sqrt(z[*,i_n])
     
  endfor

  x[*] = kr_rho[*]
  y[*] = kt_rho[*]

  plot_def_new,pname
  
  if (log_flag eq 1) then z = alog10(z)

  minz=min(z)
  maxz=max(z)

  print, 'minz=',minz,'maxz=',maxz,'nlevels=',nlevels

  print, 'kr_rho_mid=',kr_rho[n_r/2-1]

  for i_n=0,n_n-1 do begin
     print, 'kt_rho= ',kt_rho[i_n],'  z= ',z[n_r/2-1,i_n], ' z_norm= ',z[n_r/2-1,i_n]-maxz
  endfor

  print, 'kx data'

  for i=n_r/2-1, n_r-1 do begin
     print, 'kr_rho=', kr_rho[i],'  zx= ',z[i,0], ' zx_norm= ',z[i,0]-maxz
  endfor

  ;; print data

  openw,1,'xy_pow_data'
  
  printf,1, 0
  for i_n=0,n_n-1 do begin
     printf,1, kt_rho[i_n]
  endfor
  printf,1, 1
  for i_n=0,n_n-1 do begin
     printf,1, z[n_r/2-1,i_n]-maxz
  endfor
  printf,1, 2
  for i=n_r/2-1, n_r-1 do begin
     printf,1, kr_rho[i]
  endfor
  printf,1, 3
  for i=n_r/2-1, n_r-1 do begin
     printf,1, z[i,0]-maxz
  endfor

  close,1


  loadct,CTab

  contour,z,x,y, $
  nlevels=nlevels,$
  /fill, $
  xtitle='!3'+kr_rho_string, $
  xstyle=1,$
;;          xrange=[min(x),max(x)],$
  xrange=[0.0,max(y)],$
  ytitle='!3'+kt_rho_string, $
  ystyle=1,$
;;          yrange=[min(y),max(y)],$
  yrange=[0.0,max(y)],$
  title='!3'+title

  plot_finish

  set_line_colors

  return
  
end
