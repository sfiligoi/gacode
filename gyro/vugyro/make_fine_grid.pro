pro make_fine_grid

  common GLOBAL
  common POLOIDAL_DATA

  ny = n_theta_plot*mUt

  qt = fltarr(ny,n_r)
  th = fltarr(ny)

  dth = 2*!Pi/ny  
  th_min = -!Pi
  th_max = th_min+(ny-1)*dth

  for j=0,ny-1 do begin
     th[j] = th_min+dth*j
  endfor

  qt[*,*] = th # q_i

end
