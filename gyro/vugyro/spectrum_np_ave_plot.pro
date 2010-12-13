pro spectrum_np_ave_plot

  common GLOBAL
  common PLOT_VARIABLES

  ;;---------------------------------
  ;; Minimum for log-spectrum
  spec_z_min = 1e-6
  ;;---------------------------------
 
  z = fltarr(n_n+1,n_r+1)

  z[*,*] = 0.0

  for i_t=it1,it2 do begin
    z[0:n_n-1,0:n_r-1] = z[0:n_n-1,0:n_r-1]+transpose(kxkyspec[*,*,i_t])
  endfor

  if (p_squ eq 1) then begin
    for p=0,n_r-1 do begin
     z[0:n_n-1,p] = sqrt(abs(z[0:n_n-1,p]))*(p-n_r/2)^2
    endfor
    z[*,*] = z[*,*]^2
  endif

  xtitle='!3n'
  ytitle='!3p'
  title = '['+strtrim(string(max(sqrt(z))),2)+']'

  plot_def_new,'spectrum_np_ave'
  
  ;;------------------------------------------------
  ;; Construct n vector
  ;; 
  n  = n_dn*(indgen(n_n+1))+n0
  ir = (indgen(n_r+1))-n_r/2-0.5  

  if i_log eq 1 then begin

    surface,sqrt(z/max(z)),n,ir,$
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

    surface,sqrt(z/max(z)),n,ir,$
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
