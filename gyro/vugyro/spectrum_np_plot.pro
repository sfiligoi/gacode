pro spectrum_np_plot

  common GLOBAL
  common PLOT_VARIABLES

  ;;---------------------------------
  ;; Minimum for log-spectrum
  spec_z_min = 1e-6
  ;;---------------------------------
  

  z = fltarr(n_n+1,n_r+1)

  z_max = max(abs(kxkyspec[*,*,t_c]))
  z[0:n_n-1,0:n_r-1] = sqrt(transpose(kxkyspec[*,*,t_c])/z_max)
  
  if (p_squ eq 1) then begin
     for p=0,n_r-1 do begin
        z[0:n_n-1,p] = sqrt(kxkyspec[p,*,t_c])*(p-n_r/2)^2
     endfor
     z_max = max(abs(z[*,*]))
     z[*,*] = z[*,*]/z_max
  endif

  xtitle='!3n'
  ytitle='!3p'
  title = '['+strtrim(string(z_max),2)+']'

  plot_def_new,'np_phi'
  
  ;;------------------------------------------------
  ;; Construct n vector
  ;; 
  n  = n_dn*(indgen(n_n+1))+n0
  ir = (indgen(n_r+1))-n_r/2-0.5  

  if i_log eq 1 then begin

     surface,z,n,ir,$
             ax=ax_c, $
             az=az_c, $
             xtitle=xtitle, $
             ytitle=ytitle, $
             ztitle='!4u!3!dp!n!un!n',$
             title=title, $
             /lego, $
             charsize=1.4*c_size,$
             min_value=spec_z_min, $
             zrange=[0.005,1.0],$
             /zlog, $
             color=line

  endif else begin

     surface,z,n,ir,$
             ax=ax_c, $
             az=az_c, $
             xtitle=xtitle, $
             ytitle=ytitle, $
             ztitle='!4u!3!dp!n!un!n',$
             title=title, $
             /lego, $
             charsize=1.4*c_size,$
             min_value=0.0, $
             color=line

  endelse

  plot_finish

end 
