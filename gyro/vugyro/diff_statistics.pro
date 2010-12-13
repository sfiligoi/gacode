pro diff_statistics,y,t1,t2,i1,i2,A_f,D_RMS

  common GLOBAL

  t_indices,t1,t2,i1,i2

  if (i2 eq i1) then begin
    print,'Error: times too close'
    i2 = i2+1
  endif

  T_total = t[i2]-t[i1]

  sum_f  = 0.0
  sum_f2 = 0.0

  for i=i1,i2-1 do begin

    dt = t[i+1]-t[i]
    f  = 0.5*(y[i+1]+y[i])

    sum_f  = sum_f + f*dt
    sum_f2 = sum_f2 + f*f*dt
    
  endfor

  A_f  = sum_f/T_total
  A_f2 = sum_f2/T_total

  ;; Added abs() to cure case near zero
  
  D_RMS = sqrt(abs(A_f2-A_f^2))

  return

end

   

   
