pro safe_read,variable,default

  common GLOBAL

  if (not (eof(1))) then begin
    readf,1,variable
  endif else begin
    variable=default
  endelse

  return

end
