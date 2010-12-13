pro smf_tag,tag,ptag,ytag

  common GLOBAL

  if (i_spec lt n_ion) then begin
     tag_s = 'Ion '+strtrim(string(i_spec+1),2)
     ptag_s = 'ion'+strtrim(string(i_spec+1),2)
     ytag_s = '!di!n'
  endif else begin
     tag_s = 'Electron'
     ptag_s = 'elec'
     ytag_s = '!de!n'
  endelse

  if (i_moment eq 0) then begin
     tag_m = 'Density'
     ptag_m = 'density'
     ytag_m = '!3D'
  endif else begin
     tag_m = 'Energy'
     ptag_m = 'energy'
     ytag_m = '!4v!3'
  endelse

  if (i_f eq 0) then begin
     tag_f = 'ES'
     ptag_f = 'es'
  endif else if (i_f eq 1) then begin
     tag_f = 'EM'
     ptag_f = 'em'
  endif else begin
     tag_f = 'Tot'
     ptag_f = 'tot'
  endelse

  tag  = tag_f+' '+tag_s+' '+tag_m
  ptag = ptag_s+'-'+ptag_m+'-'+ptag_f
  ytag = ytag_m+ytag_s

end

