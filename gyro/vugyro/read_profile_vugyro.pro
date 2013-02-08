;-------------------------------------------------------------
; ** profile_vugyro, if not up-to-date, should be brought 
; up-to-date using "gyro -t". **
;
; Units conversion data should eventually be moved out of
; profile_vugyro.out and all conversions done from units.out.
;-------------------------------------------------------------

pro read_profile_vugyro

  common GLOBAL
  common PROFILE_SIM_DATA

  file = 'out.gyro.profile'  

  ;;------------------------------------------------
  ;; Want these set explicitly to integers
  ;;
  n_theta_section = 0
  n_fine          = 0
  n_r             = 0
  n_kinetic       = 0
  n_pass          = 0
  n_trap          = 0
  n_energy        = 0
  n_theta_plot    = 0
  n0              = 0
  n_n             = 0
  n_dn            = 0
  nonlinear_flag  = 0
  electron_method = 0
  boundary_method = 0
  ;;------------------------------------------------

  ;;------------------------------------------------
  ;; Simulation profile label definitions
  ;;
  profile_label = ['r',$
                   'q',$
                   'r_s',$
                   'q_s',$
                   'dlntdr_s',$
                   'dlnndr_s',$
                   'tem_s',$
                   'den_s',$
                   'aspect_s',$
                   'delta_s',$
                   'zeta_s',$
                   'kappa_s',$
                   'drmaj_s',$
                   'shat_s',$
                   's_delta_s',$                   
                   's_zeta_s',$                   
                   's_kappa_s',$
                   'zmag_s',$
                   'dzmag_s',$
                   'beta_unit_s',$
                   'gamma_e_s',$
                   'gamma_p_s',$
                   'mach_s',$
                   'b_unit_s',$
                   'dr_eodr',$                   
                   'z_eff_s',$                   
                   'nu_s',$
                   'w0_s', $
                   'density_imbalance']

  n_profile_label = n_elements(profile_label)
  ;;------------------------------------------------

  openr,1,file,error=i_err
  if (i_err eq 0) then begin

     ;;---------------------------------------------------
     ;; START READING
     ;;---------------------------------------------------

     ;; essential scalars
     readf,1,n_r
     readf,1,n_theta_section
     readf,1,n_pass
     readf,1,n_trap
     readf,1,n_energy
     readf,1,n_theta_plot
     readf,1,n_0
     readf,1,n_n
     readf,1,n_dn
     readf,1,n_bnd
     readf,1,nonlinear_flag
     readf,1,electron_method
     readf,1,n_field
     readf,1,n_ion
     readf,1,n_kinetic
     readf,1,n_spec
;     readf,1,field_r0_flag
;     readf,1,field_r0_grid
     readf,1,n_rho
     readf,1,boundary_method

     ;; Need to do some calculations, then allocate 
     ;; next arrays to be read

     n_lambda = n_pass+n_trap

     r           = fltarr(n_r)
     q_i         = fltarr(n_r)
     r_m         = fltarr(n_r)
     q_m         = fltarr(n_r)
     dlntdr_s    = fltarr(n_spec,n_r)
     dlnndr_s    = fltarr(n_spec,n_r)
     tem_s       = fltarr(n_spec,n_r)
     den_s       = fltarr(n_spec,n_r)
     aspect_s    = fltarr(n_r)
     delta_s     = fltarr(n_r)
     zeta_s      = fltarr(n_r)
     kappa_s     = fltarr(n_r)
     shift_s     = fltarr(n_r)
     shat_s      = fltarr(n_r)
     s_delta_s   = fltarr(n_r)
     s_zeta_s    = fltarr(n_r)
     s_kappa_s   = fltarr(n_r)
     zmag_s      = fltarr(n_r)
     dzmag_s     = fltarr(n_r)
     beta_unit_s = fltarr(n_r)
     gamma_e_s   = fltarr(n_r)
     gamma_p_s   = fltarr(n_r)
     mach_s      = fltarr(n_r)
     b_unit_s    = fltarr(n_r)
     dr_eodr     = fltarr(n_r)
     z_eff_s     = fltarr(n_r)
     nu_s        = fltarr(n_r)
     w0_s        = fltarr(n_r)

     ; don't exist anymore

     chi_i_exp       = fltarr(n_r)
     chi_e_exp       = fltarr(n_r)
     diff_to_flow_e1 = fltarr(n_r)
     diff_to_flow_e2 = fltarr(n_r)
     eta_i_tot_exp   = fltarr(n_r)
     aolvi_exp       = fltarr(n_r)
     diff_to_flow_mi = fltarr(n_r) 
     diff_ne_exp     = fltarr(n_r)
     aolne_exp       = fltarr(n_r)
     diff_to_flow_ne = fltarr(n_r)
     eta_i_diff_exp  = fltarr(n_r)
     diff_to_flow_heating = fltarr(n_r)

     lambda  = fltarr(n_lambda)
     energy  = fltarr(n_energy)
     kt_rho  = fltarr(n_n)
     zcharge = fltarr(n_spec)

     readf,1,r
     readf,1,q_i
     readf,1,r_m
     readf,1,q_m
     readf,1,dlntdr_s
     readf,1,dlnndr_s
     readf,1,tem_s
     readf,1,den_s
     readf,1,aspect_s
     readf,1,delta_s
     readf,1,zeta_s
     readf,1,kappa_s
     readf,1,shift_s
     readf,1,shat_s
     readf,1,s_delta_s
     readf,1,s_zeta_s 
     readf,1,s_kappa_s
     readf,1,zmag_s
     readf,1,dzmag_s
     readf,1,beta_unit_s
     readf,1,gamma_e_s
     readf,1,gamma_p_s
     readf,1,mach_s
     readf,1,b_unit_s
     readf,1,dr_eodr
     readf,1,z_eff_s
     readf,1,nu_s
     readf,1,w0_s

     readf,1,box_multiplier

     ;; velocity space, etc
     readf,1,lambda
     readf,1,energy
     readf,1,lambda_tp     
     readf,1,kt_rho
     readf,1,rho_s
     readf,1,zcharge

     ;; n_fine version control
     readf,1,n_fine 

     ;;---------------------------------------------------
     ;; FINISHED READING
     ;;---------------------------------------------------

     for ii=0,n_r-1 do begin
        eta_i_diff_exp[ii] = eta_i_tot_exp[ii]- $
          diff_ne_exp[ii]*aolne_exp[ii]/aolvi_exp[ii]
     endfor
     
     ;; Basic grid initializations

     n_tor = n_0+indgen(n_n)*n_dn

     ;; Initialize singular surface counter

     n_ss = 0
     if (n_n gt 1) then begin
        dn_ss = n_tor[1]
     endif else begin 
        dn_ss = n_tor[0]
     endelse

     ;; Initialize plotting vectors

     if n_theta_plot eq 1 then begin

        theta_plot = [0.0]
        thetad_plot = [0.0]

     endif else begin

        theta_plot = -!pi+2*!pi*findgen(n_theta_plot)/n_theta_plot
        thetad_plot = -1.0+2*findgen(n_theta_plot)/n_theta_plot

     endelse

     ;; Initialize plotting vectors

;     if field_r0_grid eq 1 then begin

;        thetad_r0_plot = [0.0]

;     endif else begin

;        thetad_r0_plot = -1.0+2*findgen(field_r0_grid)/field_r0_grid

;     endelse

     ;; Initialize counter for poloidal harmonics

     m_c = q_i[n_r/2]*n_tor[n_n/2]

     ;; Initialize distribution function dimensions

     n_stack   = 4*n_theta_section-4
     n_theta_p = 2*n_theta_section-2
     n_theta_t = 2*n_theta_section-1

     j_c = n_theta_plot/2

     spec_label = strarr(n_spec)
     for i=0,n_spec-2 do begin
        spec_label[i] = 'ion_'+strtrim(string(i+1),2)+'_'
     endfor

     spec_label[n_spec-1] = 'electron_'

     print_found,file,1
     exists_profile = 1

  endif else begin

     print,'FATAL ERROR: missing profile_vugyro.out'
     stop

  endelse

  close,1

  return

end  
