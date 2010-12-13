pro get_singsurf_vec,n,r_surf

  common GLOBAL

  root = fltarr(10*n_r)

  if (n eq 0) then begin
     i_root = 0
     goto, finish
  endif

  q_diff = fltarr(n_r)

  i_root = 0
  for m=1,10*n do begin

     q_0 = float(m)/n

     q_diff[*] = q_i[*]-q_0

     for i=0,n_r-2 do begin
        if (q_diff[i] le 0.0 and q_diff[i+1] ge 0.0) then begin
           root[i_root] = $
             r[i]+(0.0-q_diff[i])/(q_diff[i+1]-q_diff[i])*(r[i+1]-r[i])
           i_root = i_root+1
        endif 
        if (q_diff[i] ge 0.0 and q_diff[i+1] le 0.0) then begin
           root[i_root] = $
             r[i]+(0.0-q_diff[i])/(q_diff[i+1]-q_diff[i])*(r[i+1]-r[i])
           i_root = i_root+1
        endif
     endfor 
  endfor

  finish:

  if i_root eq 0 then begin 
     i_root = 1
     root[0] = 0.0
  endif

  r_surf = fltarr(i_root)
  for i=0,i_root-1 do begin
     r_surf[i] = root[i]
  endfor

end
