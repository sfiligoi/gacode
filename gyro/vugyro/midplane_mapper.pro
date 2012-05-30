pro midplane_mapper
;
; 09 Nov 2007: CH
;  added in temperature flucuations if both density and energy moments
;  available, define T = (2/3)E - n
;
  common GLOBAL
  common MIDPLANE_DATA

  n_pwr = n_field*exists_u+(exists_mom_n+exists_mom_e+exists_mom_v)*n_kinetic

  if (exists_mom_n AND exists_mom_e) then begin
      n_pwr += n_kinetic
      print, 'adding temp flucs to midplane power'
  endif

  rp = indgen(n_r)-(n_r/2-1)
  kr_rho = rp*2*!pi/(r[n_r-1]-r[0])*rho_s

  if n_pwr eq 0 then return

; release previous allocation of pwr:
  pwr=0
  pwr = complexarr(n_pwr,n_r,n_n,n_time)

  tag_spec = strarr(n_kinetic)
  tag      = strarr(n_pwr)

  for ispec=0, n_kinetic-1 do tag_spec[ispec] = 'ion'+$
                               string(ispec+1,format='(I1)')
; correct the last label if there are electrons:
  if (n_kinetic gt 1) then tag_spec[n_kinetic-1] = 'elec'

  if n_theta_plot eq 1 then begin
     j = 0
  endif else begin 
     j = n_theta_plot/2
  endelse

  i_pwr = 0
  if exists_u eq 1 then begin
     pwr[i_pwr,*,*,*] = complex(u[0,j,*,0,*,*],$
                                u[1,j,*,0,*,*])
     tag[i_pwr] = 'potential'
     i_pwr = i_pwr+1
     if n_field eq 2 then begin
     pwr[i_pwr,*,*,*] = complex(u[0,j,*,1,*,*],$
                                u[1,j,*,1,*,*])
     tag[i_pwr] = 'A_par'
     i_pwr = i_pwr+1
     endif
  endif

  if exists_mom_n eq 1 then begin
     for is=0,n_kinetic-1 do begin
        pwr[i_pwr,*,*,*] = complex(mom_n[0,j,*,is,*,*],$
                                   mom_n[1,j,*,is,*,*])
        tag[i_pwr] = 'density_'+tag_spec[is]
        i_pwr = i_pwr+1
     endfor
  endif

  if exists_mom_e eq 1 then begin
     for is=0,n_kinetic-1 do begin
        pwr[i_pwr,*,*,*] = complex(mom_e[0,j,*,is,*,*],$
                                   mom_e[1,j,*,is,*,*])
        tag[i_pwr] = 'energy_'+tag_spec[is]
        i_pwr = i_pwr+1
     endfor
  endif

  if exists_mom_v eq 1 then begin
     for is=0,n_kinetic-1 do begin
        pwr[i_pwr,*,*,*] = complex(mom_v[0,j,*,is,*,*],$
                                   mom_v[1,j,*,is,*,*])
        tag[i_pwr] = 'vparallel_'+tag_spec[is]
        i_pwr = i_pwr+1
     endfor
  endif

  if ((exists_mom_e eq 1) AND (exists_mom_n eq 1)) then begin
     for is=0,n_kinetic-1 do begin
;  define E = (3/2)nT -> dE = (3/2)*(dn*T0 + n0*dT)
;  -> dT = [(2/3)dE -T0*dn]/n0
         mid_e = REFORM(complex(mom_e[0,j,*,is,*,*],$
                                mom_e[1,j,*,is,*,*]))
         mid_n = REFORM(complex(mom_n[0,j,*,is,*,*],$
                                mom_n[1,j,*,is,*,*]))

         for ir = 0, n_r-1 do pwr[i_pwr,ir,*,*] = $
           REFORM(((2./3)*mid_e[ir,*,*] - tem_s[is,ir]*mid_n[ir,*,*])/$
                  den_s[is,ir])

         tag[i_pwr] = 'tem_'+tag_spec[is]
         i_pwr = i_pwr+1
     endfor
  endif

  i_pwr = 0

end



