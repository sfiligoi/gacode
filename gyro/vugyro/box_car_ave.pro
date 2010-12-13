pro box_car_ave,y,i1,i2,i_del,s_A,s_D_RMS

  common GLOBAL

  if (i2 eq i1) then begin
    print,'Error: times too close'
    i2 = i2+1
  endif


  sum_f  = 0.0
  for i=i1,i1+i_del do begin
    f  = y[i]
    sum_f  = sum_f + f
  endfor
  
  y_bca = fltarr(n_time)
  y_bca[*] = 0.
  for i=i1,i1+i_del do begin
   y_bca(i) = sum_f
  endfor

  for i=i1+i_del,i2-1 do begin
   y_bca(i+1)=y_bca(i)+(y(i+1)-y(i-i_del))
  endfor

  y_bca[*]=y_bca[*]/i_del

  y[*]= y_bca[*]

 ;;get average and D_RMS


   i_total = i2-(i1+i_del)

  sum_f  = 0.0
  sum_f2 = 0.0

  for i=i1+i_del,i2-1 do begin

    f  = 0.5*(y_bca[i+1]+y_bca[i])

    sum_f = sum_f + f
    sum_f2 = sum_f2 + f*f

  endfor

  A_f  = sum_f/i_total
  A_f2 = sum_f2/i_total

  D_RMS = sqrt(abs(A_f2-A_f^2))

  print, '  '
  print, ' bca_A=', A_f,  ' bca_D_RMS=', D_RMS
  

 ;; samples

  i=i1+i_del
  i_cnt_max=(i2-i1)/i_del
   print, '   '
   sum_xf=0.
   sum_xf2=0.
   print, '  '
   print, ' i1=',i1,' i2=',i2,' n_time=',n_time
   print, '  '
  for i_cnt = 1, i_cnt_max do begin
  ;; if (i ge 0 and i lt n_time-1) then begin
    if (i ge 0 and i lt n_time) then begin
    print, ' i_cnt=', i_cnt, ' i=',i, ' y_bca=',y_bca(i)
    xf=y_bca(i)
    sum_xf = sum_xf + xf
    sum_xf2 = sum_xf2 + xf*xf
    i_cnt_last = i_cnt
   endif
     i=i+i_del
  endfor

  s_A  = sum_xf/i_cnt_last
  s_A_2 = sum_xf2/i_cnt_last
  s_D_RMS = sqrt(abs(s_A_2-s_A^2))

  print, '   '
  print, ' s_A=', s_A,  ' s_D_RMS=', s_D_RMS
  
  return

  end
   
