pro tag_gen,n_ion,i_spec,i_moment,i_f,i_profile

  common TAG_DATA

  if (i_spec lt n_ion) then begin
     tag_s = 'Ion '+strtrim(string(i_spec+1),2)
     ftag_s = 'ion'+strtrim(string(i_spec+1),2)
     ltag_s = '!di!n'
  endif else begin
     tag_s = 'Electron'
     ftag_s = 'elec'
     ltag_s = '!de!n'
  endelse

  if (i_moment eq 0) then begin
     if i_profile eq 0 then begin
        tag_m = 'Density'
        ftag_m = 'density'
        ltag_m = '!3D'
     endif else begin 
        tag_m = 'Density'
        ftag_m = 'density'
        ltag_m = '!3n'
     endelse
  endif else begin
     if i_profile eq 0 then begin
        tag_m = 'Energy'
        ftag_m = 'energy'
        ltag_m = '!4v!3'
     endif else begin 
        tag_m = 'Temperature'
        ftag_m = 'temperature'
        ltag_m = '!3T'
     endelse
  endelse

  if (i_f eq 0) then begin
     tag_f = 'ES'
     ftag_f = 'es'
  endif else if (i_f eq 1) then begin
     tag_f = 'EM'
     ftag_f = 'em'
  endif else begin
     tag_f = 'Tot'
     ftag_f = 'tot'
  endelse

end

