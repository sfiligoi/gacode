pro get_singsurf_den,r_surf

  common GLOBAL
  common PRIVATE_OSC

  ;; singular surface density calculation

  r_all_surf = fltarr(100000)
  n_cnt=0
  for n_level=1,n_n-1 do begin
     get_singsurf_vec,n_level*dn_ss,r_surf
     n_surf = n_elements(r_surf)

     for i=0,n_surf-1 do begin
        r_all_surf(i+n_cnt)=r_surf(i)
     endfor
     n_cnt = n_cnt+n_surf
  endfor

  print , 'n_cnt=',n_cnt

  dup = fltarr(n_cnt)
  dup[*]=0.
  for i=0,n_cnt-1 do begin
     for ii = 0,n_cnt-1 do begin
        if r_all_surf(i) eq r_all_surf(ii) then begin
           dup(i) = dup(i)+1.
        endif
     endfor
  endfor

  tot_surf = 0.
  for i=0,n_cnt-1 do begin
     tot_surf=tot_surf+1./dup(i)
  endfor

  print , '  '
  print , 'tot_surf=',tot_surf 

  del_ss_bin = (r(n_r-1)-r(0))/n_ss_bin
  den_ss = fltarr(n_ss_bin)
  den_sn = fltarr(n_ss_bin)

  den_ss[*]=0.
  den_sn[*]=0.
  tot_ss = 0.
  tot_mn = 0.
  tot_ss_den = 0.
  tot_mn_den = 0.
  for n_bin = 0,n_ss_bin-1 do begin
     for i=0,n_cnt-1 do begin
        if n_bin*del_ss_bin lt r_all_surf(i)-r(0) then begin
           if (n_bin+1)*del_ss_bin ge r_all_surf(i)-r(0) then begin
              den_ss(n_bin) = den_ss(n_bin) +1./dup(i)
              den_sn(n_bin) = den_sn(n_bin) + 1
              tot_ss = tot_ss + 1./dup(i)
              tot_mn = tot_mn +1.
           endif
        endif
     endfor
     tot_ss_den = tot_ss_den + den_ss(n_bin)
     tot_mn_den = tot_mn_den + den_sn(n_bin)
  endfor
  print , 'tot_ss=',tot_ss,' tot_mn=',tot_mn
  print , 'tot_ss_den=',tot_ss_den, ' tot_mn_den=',tot_mn_den
  for n_bin = 0,n_ss_bin-1 do begin
     den_ss(n_bin)=n_ss_bin*den_ss(n_bin)/tot_ss_den
     den_sn(n_bin)=n_ss_bin*den_sn(n_bin)/tot_mn_den
  endfor

  plot_den_ss=fltarr(2*n_bin)
  plot_den_mn=fltarr(2*n_bin)
  plot_r_den_ss=fltarr(2*n_bin)
  nn_bin = 0
  for n_bin = 0,n_ss_bin-1 do begin
     nn_bin = 2*n_bin
     plot_den_ss(nn_bin)=den_ss(n_bin)
     plot_den_mn(nn_bin)=den_sn(n_bin)
     plot_r_den_ss(nn_bin)=r(0)+n_bin*del_ss_bin
     nn_bin=nn_bin+1
     plot_den_ss(nn_bin)=den_ss(n_bin)
     plot_den_mn(nn_bin)=den_sn(n_bin)
     plot_r_den_ss(nn_bin)=r(0)+(n_bin+1)*del_ss_bin-0.0001
  endfor

  ;; end singular surface density calulation

  if i_ss_mn eq 0 then begin
     oplot,plot_r_den_ss,plot_den_ss,color=line
     print, 'density of singular surfaces i_ss_mn = 0'
  end else begin
     oplot,plot_r_den_ss,plot_den_mn,color=line
     print, 'density of mn modes i_ss_mn = 1'
  endelse

  ;; end singular surface density and mn mode density plot
  
end
