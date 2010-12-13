pro print_found,file,flag

  if flag eq 1 then begin
     print,'* READ: '+file
  endif else begin
     print,'* READ: '+file+' [MISSING]'
  endelse

end
