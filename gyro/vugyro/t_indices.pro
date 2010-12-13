pro t_indices,t1,t2,i1,i2

  common GLOBAL

  if (t1 ge t2) or $
      (t1 lt 0)  or $
      (t2 lt 0)  or $
      (t1 gt t[n_time1]) or $
      (t2 gt t[n_time1]) then begin

    print,'Correcting time-average interval.'

    t1 = t[0]
    t2 = t[n_time1]

  endif

  i = 0
  while (t[i] lt t1) do i=i+1
  i1 = i
  while (t[i] lt t2) do i=i+1
  i2 = i

  return

end

   

   
