pro vector_deriv,x,f,fx,bc

  ; periodic    : bc = 1
  ; nonperiodic : bc = 2

  n = n_elements(x)

  if (bc eq 2) then begin

     ; nonperiodic

     fx[0,*] = (f[1,*]-f[0,*])/(x[1]-x[0])
     for i=1,n-2 do begin
        fx[i,*] = (f[i+1,*]-f[i-1,*])/(x[i+1]-x[i-1])
     endfor
     fx[n-1,*] = (f[n-1,*]-f[n-2,*])/(x[n-1]-x[n-2])

  endif else begin

     ; periodic

     dx = x[2]-x[0]

     fx[0,*] = (f[1,*]-f[n-1,*])/dx
     for i=1,n-2 do begin
        fx[i,*] = (f[i+1,*]-f[i-1,*])/dx
     endfor
     fx[n-1,*] = (f[0,*]-f[n-2,*])/dx

  endelse

  return

end
