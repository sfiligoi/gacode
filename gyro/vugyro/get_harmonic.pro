pro get_harmonic,f_m,r_surf

  common GLOBAL

  ;; project (m_c)-poloidal harmonic 

  f_m = complexarr(n_r)
  r_surf = fltarr(2)

  dtheta  = 2*!pi/n_theta_plot

  phi_j = complexarr(n_theta_plot,n_r)

  phi_j[*,*] = complex(u[0,*,*,i_f,in_c,t_c],u[1,*,*,i_f,in_c,t_c])

  s  = complexarr(n_r)
  z  = complexarr(n_r)
  w  = complexarr(n_r)
  ex = complexarr(n_r)
  
  s  = complex(0,1)*(n_tor[in_c]*q_i-m_c)
  z  = s*dtheta
  ex = exp(z)

  for i=0,n_r-1 do begin

     if abs(z[i]) lt 1e-2 then begin
        
        w[i] = 1.0

     endif else begin

        w[i] = (ex[i]-1.0)/z[i]+(1/ex[i]-1.0)*(1.0+(z[i]-1.0)*ex[i])/z[i]^2

     endelse

     sum = 0.0
     for m=0,n_theta_plot-1 do begin
        sum = sum + exp(s[i]*theta_plot[m])*phi_j[m,i]
     endfor

     f_m[i] = sum*w[i]/n_theta_plot

  endfor

  ;; Find singular surface(s) (r_surf[0],r_surf[1])

  q_center = float(m_c)/float(n_tor[in_c])

  q_diff = q_i-q_center

  for i=0,n_r-2 do begin
     if (q_diff[i] le 0.0 and q_diff[i+1] ge 0.0) then begin
        r_surf[0] = $
          r[i]+(0.0-q_diff[i])/(q_diff[i+1]-q_diff[i])*(r[i+1]-r[i])
     endif
     if (q_diff[i] ge 0.0 and q_diff[i+1] le 0.0) then begin
        r_surf[1] = $
          r[i]+(0.0-q_diff[i])/(q_diff[i+1]-q_diff[i])*(r[i+1]-r[i])
     endif
  endfor 

  print,r_surf

end
