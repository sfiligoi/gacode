; The coding below reconstructs a poloidal 
; map of the real potential

pro extract_contour,tt

  common GLOBAL
  common POLOIDAL_DATA

  ny = n_theta_plot*mUt
  a = fltarr(ny,n_r)

  openw, 1, 'r_theta_nu.txt'
  for irad=0,n_r-1 do begin
     for ith=0,n_theta_plot-1 do begin
       nu=geometry[0,ith,irad]
       printf,1,'{',r[irad],',',th[ith],',',nu,'}'
     endfor
  endfor
  close,1

  for nn=0,n_n-1 do begin
    nst = strtrim(string(fix(n_tor[nn])),2)
    tst = strtrim(string(fix(t[tt])),2)
    pname = 'phi_Apar_' + nst + '_' + tst + '.txt'
    openw, 1, pname
    for irad=0,n_r-1 do begin
       for ith=0,n_theta_plot-1 do begin
;         nu=geometry[0,ith,irad]
         phir=U[0,ith,irad,0,nn,tt]
         phii=U[1,ith,irad,0,nn,tt]
         Apr =U[0,ith,irad,1,nn,tt]
         Api =U[1,ith,irad,1,nn,tt]
         printf,1,'{',phir,',',phii,',',Apr,',',Api,'}'
       endfor
    endfor
    close,1
  endfor


  if i_field le 4 then begin

    U_temp = complexarr(n_theta_plot+1,n_r)
    U_c    = complexarr(ny,n_r,n_n)
    U_int  = complexarr(ny+1)
    phase  = complexarr(n_r)

    for nn=0,n_n-1 do begin

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

        4: begin

          U_temp[0:n_theta_plot-1,*] = $
              complex(mom_v[0,0:n_theta_plot-1,*,i_spec,nn,tt], $
                      mom_v[1,0:n_theta_plot-1,*,i_spec,nn,tt])

          U_temp[n_theta_plot,*] = $
              complex(mom_v[0,0,*,i_spec,nn,tt], $
                      mom_v[1,0,*,i_spec,nn,tt])*phase[*]

        end

      endcase

      for i=0,n_r-1 do begin
        U_int = interpol(U_temp[*,i],ny+1)
        U_c[0:ny-1,i,nn] = U_int[0:ny-1]
      endfor

    endfor
    
    if (n_tor[0] eq 0) then begin

      ;; With zero mode;

      if ave_poloidal eq 0 then begin

        a[*,*] = float(U_c[*,*,0])

      endif else begin

        for i=0,n_r-1 do begin
          U0_ave = total(float(U_c[*,i,0]))/ny
          a[*,i] = float(U_c[*,i,0])-U0_ave
        endfor

      endelse

      for nn=1,n_n-1 do begin
        a[*,*] = a[*,*] + 2*float(U_c[*,*,nn]* $
                                  exp(complex(0,1)*n_tor[nn]*qt[*,*]))
      endfor

    endif else begin

      ;; Without zero mode;

      a[*,*] = 0.0
      for nn=0,n_n-1 do begin
        a[*,*] = a[*,*] + 2*float(U_c[*,*,nn]* $
                                  exp(complex(0,1)*n_tor[nn]*qt[*,*]))
      endfor

    endelse

  endif else begin              

    ;;**************************************;

    a1 = fltarr(ny,n_r)
    a2 = fltarr(ny,n_r)

    U1_real = fltarr(ny,n_r)
    U1_imag = fltarr(ny,n_r)
    U1_temp = fltarr(n_theta_plot,n_r)
    U1_c    = complexarr(ny,n_r,n_n)
    U2_real = fltarr(ny,n_r)
    U2_imag = fltarr(ny,n_r)
    U2_temp = fltarr(n_theta_plot,n_r)
    U2_c    = complexarr(ny,n_r,n_n)

    case (i_field) of

      4: begin

        for nn=0,n_n-1 do begin
          U1_temp[*,*] = U[0,*,*,0,nn,tt]
          U1_real = rebin(U1_temp,ny,n_r)
          U1_temp[*,*] = U[1,*,*,0,nn,tt]
          U1_imag = rebin(U1_temp,ny,n_r)
          U1_c[*,*,nn] = complex(U1_real[*,*],U1_imag[*,*])
        endfor

        for nn=0,n_n-1 do begin
          U2_temp[*,*] = mom_n[0,*,*,i_spec,nn,tt]
          U2_real = rebin(U2_temp,ny,n_r)
          U2_temp[*,*] = mom_n[1,*,*,i_spec,nn,tt]
          U2_imag = rebin(U2_temp,ny,n_r)
          U2_c[*,*,nn] = complex(U2_real[*,*],U2_imag[*,*])
        endfor

      end

      5: begin

        for nn=0,n_n-1 do begin
          U1_temp[*,*] = U[0,*,*,0,nn,tt]
          U1_real = rebin(U1_temp,ny,n_r)
          U1_temp[*,*] = U[1,*,*,0,nn,tt]
          U1_imag = rebin(U1_temp,ny,n_r)
          U1_c[*,*,nn] = complex(U1_real[*,*],U1_imag[*,*])
        endfor

        for nn=0,n_n-1 do begin
          U2_temp[*,*] = mom_e[0,*,*,i_spec,nn,tt]
          U2_real = rebin(U2_temp,ny,n_r)
          U2_temp[*,*] = mom_e[1,*,*,i_spec,nn,tt]
          U2_imag = rebin(U2_temp,ny,n_r)
          U2_c[*,*,nn] = complex(U2_real[*,*],U2_imag[*,*])
        endfor

      end

    endcase    

    ;; Without zero mode;

    a1[*,*] = 0.0
    for nn=0,n_n-1 do begin
      a1[*,*] = a1[*,*] + 2*float(complex(0,1)*n_tor[nn]*U1_c[*,*,nn]* $
                                  exp(complex(0,1)*n_tor[nn]*qt[*,*]))
    endfor
    a2[*,*] = 0.0
    for nn=0,n_n-1 do begin
      a2[*,*] = a2[*,*] + 2*float(U2_c[*,*,nn]* $
                                  exp(complex(0,1)*n_tor[nn]*qt[*,*]))
    endfor

    a = -a1*a2

  endelse   

end
