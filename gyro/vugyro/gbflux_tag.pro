pro gbflux_tag,tag,ptag,ytag,units1,units2,units3

  common GLOBAL

  if (i_spec lt n_ion) and (electron_method ne 3) then begin
     tag_s  = 'ion'+strtrim(string(i_spec+1),2)
     ptag_s = 'ion'+strtrim(string(i_spec+1),2)
     ytag_s = '!di!n'
  endif else begin
     tag_s  = 'elec'
     ptag_s = 'elec'
     ytag_s = '!de!n'
  endelse

  case i_moment of

     0: begin  
        tag_m  = 'density'
        ptag_m = 'density'
        ytag_m = '!4C!3'
        units1  = '/!4C!3!dGB!n'
        units2  = ' MW/keV/m!u2!n'
        units3  = '(dV/dr) MW/keV'
    end

     1: begin
        tag_m  = 'energy'
        ptag_m = 'energy'
        ytag_m = 'Q'
        units1  = '/Q!dGB!n'
        units2  = ' MW/m!u2!n'
        units3  = '(dV/dr) MW'
     end
     
     2: begin
        tag_m  = 'momentum'
        ptag_m = 'momentum'
        ytag_m = '!4P!3'
        units1  = '/!4P!3!dGB!n'
        units2  = ' Nm/m!u2!n'
        units3  = '(dV/dr) Nm'
    end

     3: begin
        tag_m  = 'exchange'
        ptag_m = 'exchange'
        ytag_m = 'S'
        units1  = '/S!dGB!n'
        units2  = ' MW/m!u3!n'
        units3  = '(dV/dr) MW/m'
     end
     
  endcase

  case i_f of

     -1: begin
        tag_f  = 'Total'
        ptag_f = 'tot'
     end

     0: begin
        tag_f  = '!4u!3'
        ptag_f = 'phi'
     end

     1: begin
        tag_f  = 'A!d!9#!3!n'
        ptag_f = 'apar'
     end

     2: begin
        tag_f  = 'B!d!9#!3!n'
        ptag_f = 'bpar'
     end

  endcase

  tag  = '('+tag_f+','+tag_s+','+tag_m+')'
  ptag = ptag_s+'-'+ptag_m+'-'+ptag_f
  ytag = ytag_m+ytag_s

end

