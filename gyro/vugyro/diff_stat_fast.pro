pro diff_stat_fast,y,i1,i2,ave_y

  common GLOBAL

  t_total = t[i2]-t[i1]

  sum_y = 0.0

  for i=i1,i2-1 do begin

    dt = t[i+1]-t[i]
    f  = 0.5*(y[i+1]+y[i])

    sum_y = sum_y + f*dt

  endfor

  ave_y = sum_y/t_total
  
  return

end

   

   
