pro get_singsurf,r_surf

  common GLOBAL

  ;; project (m_c)-poloidal harmonic 

  r_surf = fltarr(2)

  ;; Find singular surface(s) (r_surf[0],r_surf[1])

  if n_tor[in_c] ne 0 then begin 
     q_center = m_c/float(n_tor[in_c])
  endif else begin
     q_center = -1.0
  endelse

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

  print,'r_surf=',r_surf,$
        'q_surf=',q_center,$
        'm=',m_c,$
        'n=',n_tor[in_c],$
        'in_c=',in_c
 ;;       'mf=',float(m_c),$
 ;;       'nf=',float(n_tor[in_c])

end
