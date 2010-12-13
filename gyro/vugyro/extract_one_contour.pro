; The coding below reconstructs a poloidal 
; map of the real potential for a single n

pro extract_one_contour,nn,tt

  common GLOBAL
  common POLOIDAL_DATA

  ny = n_theta_plot*mUt
  a = fltarr(ny,n_r)

  U_temp = complexarr(n_theta_plot+1,n_r)
  U_c    = complexarr(ny,n_r)
  U_int  = complexarr(ny+1)
  phase  = complexarr(n_r)

  phase[*] = exp(-2*!pi*complex(0,1)*n_tor[nn]*q_i[*])

  case (i_field) of

    1: begin

      U_temp[0:n_theta_plot-1,*] = $
          complex(U[0,0:n_theta_plot-1,*,i_f,nn,tt], $
                  U[1,0:n_theta_plot-1,*,i_f,nn,tt])

      U_temp[n_theta_plot,*] = $
          complex(U[0,0,*,i_f,nn,tt], $
                  U[1,0,*,i_f,nn,tt])*phase[*]

    end

    2: begin

      U_temp[0:n_theta_plot-1,*] = $
          complex(mom_n[0,0:n_theta_plot-1,*,i_spec,nn,tt], $
                  mom_n[1,0:n_theta_plot-1,*,i_spec,nn,tt])

      U_temp[n_theta_plot,*] = $
          complex(mom_n[0,0,*,i_spec,nn,tt], $
                  mom_n[1,0,*,i_spec,nn,tt])*phase[*]

    end

    3: begin

      U_temp[0:n_theta_plot-1,*] = $
          complex(mom_e[0,0:n_theta_plot-1,*,i_spec,nn,tt], $
                  mom_e[1,0:n_theta_plot-1,*,i_spec,nn,tt])

      U_temp[n_theta_plot,*] = $
          complex(mom_e[0,0,*,i_spec,nn,tt], $
                  mom_e[1,0,*,i_spec,nn,tt])*phase[*]

    end

  endcase

  for i=0,n_r-1 do begin
    U_int = interpol(U_temp[*,i],ny+1)
    U_c[0:ny-1,i] = U_int[0:ny-1]
  endfor

  if (n_tor[nn] eq 0) then begin

    a[*,*] = float(U_c[*,*])

  endif else begin

    a[*,*] = 2*float(U_c[*,*]* $
                     exp(complex(0,1)*n_tor[nn]*qt[*,*]))
  endelse

end
