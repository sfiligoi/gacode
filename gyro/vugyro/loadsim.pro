pro loadsim, input_dir

  common GLOBAL
  common ZMOMENT_DATA
  common T_ERROR_DATA

  widget_control,/hourglass
  
  simdir = input_dir
  print,'Changed active directory to '+simdir
  
  cd,simdir

  ;; Version tag information
  
  version_tag  = '(version not found)'
  version_date = '(time not found)' 
  openr,1,simdir+'/VERSION_tag',error=i_err

  if (i_err eq 0) then begin
     readf,1,version_tag
     readf,1,version_date
  endif else begin
     print,'No VERSION_tag found.  Setting to 1.0.0'
     version_tag = '1.0.0'
  endelse 
  close,1

  ;;-------------------------------------------------------
  ;; Read time vector;
  
  read_time_vector
  
  if (n_time gt 0) then begin
     print,'NUMBER OF TIME POINTS: ',strtrim(string(n_time),2)
  endif
  ;;-------------------------------------------------------
  
  ;;-------------------------------------------------------
  ;; RADIAL PROFILES:
  
  read_profile_vugyro
  
  exists_exp_profile = 0
  exists_exp_derived = 0

  ;; This check is needed in case a non-profile run 
  ;; is being executed in a directory which contains
  ;; profile data.

  if (n_rho gt 0) then begin

     ;; Read list of input experimental profile

     openr,1,'input.profiles.gen',err=i_err

     if (i_err eq 0) then begin

        print_found,'input.profiles.gen',1
        exists_exp_profile = 1

        readf,1,ncol_exp
        readf,1,nblock_exp
        readf,1,n_rho
        readf,1,bt_exp 
        readf,1,arho_exp

        n_exp_profile = ncol_exp*nblock_exp
        exp_profile = fltarr(n_rho,n_exp_profile)
        readf,1,exp_profile
        map_exp_profiles

     endif

     close,1

     ;; Read list of derived experimental profile

     n_exp_derived = 25

     exp_derived = fltarr(n_rho,n_exp_derived)
     read_array,exp_derived,'input.profiles.extra',exists_exp_derived

     ;; For consistency with exp_profile, use the transpose
     exp_derived = transpose(exp_derived)

  endif

  ;; Read GYRO simulation units data
  read_units

  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Read arrays of Miller geometry coefficients:
  ;;
  geometry_label = ['GEO_nu',$
                    'GEO_gsin',$
                    'GEO_gcos1',$
                    'GEO_gcos2',$
                    'GEO_usin',$
                    'GEO_ucos',$
                    'GEO_b',$
                    'GEO_g_theta',$
                    'GEO_grad_r',$
                    'GEO_g_q',$
                    'GEO_captheta']

  n_geometry = n_elements(geometry_label)

  geometry = fltarr(n_geometry,n_fine,n_r)
  read_array,geometry,'geometry_arrays.out',exists_geometry
  ;;-------------------------------------------------------;

  ;;-------------------------------------------------------
  ;; Collision-specific data
  read_collision
  ;;-------------------------------------------------------

  if n_time eq 0 then goto, no_time

  if exists_exp_profile eq 1 then begin

     ;;--------------------------------------------
     ;; Transport power flows and profiles
     ;;
     ;; 0 -> r_p(1:n_trho)       midplane minor radius 0->1
     ;; 1 -> powp_sim(1:n_trho)  simulated particle flow power (MW/Kev)
     ;; 2 -> powe_sim(1:n_trho)  simulated electron power (MW)
     ;; 3 -> powi_sim(1:n_trho)  simulated ion power (MW)
     ;; 4 -> ne_tr(1:n_trho)     transported electron density 10**19 m**(-3)
     ;; 5 -> te_tr(1:n_trho)     transported electron temperature (Kev)
     ;; 6 -> ti_tr(1:n_trho)     transported ion temperature (Kev)

     n_trho = n_rho-1
     transport_ne_te_ti= fltarr(7,n_trho,n_time)
     read_array,transport_ne_te_ti,$
                'transport_ne_te_ti.out',exists_transport_ne_te_ti
     ;;--------------------------------------------

  endif

  ;;-------------------------------------------------------
  ;; Read time-integration error 
  ;;
  t_error = fltarr(n_kinetic,n_time)
  read_array,t_error,'error.out',exists_t_error
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Read entropy
  ;;
  entropy = fltarr(n_kinetic,5,n_time)
  read_array,entropy,'entropy.out',exists_entropy
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Read velocity-space 
  
  f_nek = fltarr(n_energy,n_lambda,n_kinetic,n_field,2,n_n,n_time)
  read_array,f_nek,'flux_velocity.out',exists_velocity
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Read (kx,ky)=(kr,n) spectral data
  
  kxkyspec = fltarr(n_r,n_n,n_time)
  read_array,kxkyspec,'kxkyspec.out',exists_kxkyspec
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Read fields at r=r0:
  ;;
  field_r0 = fltarr(2,field_r0_grid,n_field,n_n,n_time)
  read_array,field_r0,'field_r0.out',exists_field_r0 
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Read large files (fields):
  ;;
  if skip_large eq 0 then begin

     u = fltarr(2,n_theta_plot,n_r,n_field,n_n,n_time)
     read_array,u,'u.out',exists_u
     
     mom_n = fltarr(2,n_theta_plot,n_r,n_kinetic,n_n,n_time)
     read_array,mom_n,'moment_n.out',exists_mom_n

     mom_e = fltarr(2,n_theta_plot,n_r,n_kinetic,n_n,n_time)
     read_array,mom_e,'moment_e.out',exists_mom_e

     mom_v = fltarr(2,n_theta_plot,n_r,n_kinetic,n_n,n_time)
     read_array,mom_v,'moment_v.out',exists_mom_v

  endif else if skip_large eq 2 then begin

     u = fltarr(2,n_theta_plot,n_r,n_field,n_n,n_time)
     read_array,u,'u.out',exists_u
     exists_mom_n = 0
     exists_mom_e = 0
     exists_mom_v = 0

  endif else begin

     exists_u = 0
     exists_mom_n = 0
     exists_mom_e = 0
     exists_mom_v = 0

  endelse
  ;;-------------------------------------------------------

  ;;-------------------------------------------------------
  ;; Frequency trace read (for linear sim only)
  ;;
  ;; 0 -> Re(Omega)
  ;; 1 -> Re(dOmega)
  ;; 2 -> Im(Omega)
  ;; 3 -> Im(dOmega)
  
  w = fltarr(4,n_n,n_time)
  read_array,w,'freq.out',exists_omega
  ;;------------------------------------------------------

  ;;------------------------------------------------------
  ;; Diffusivities:

  n_moment = 2

  diff = fltarr(n_kinetic,n_field,n_moment,n_time)
  read_array,diff,'diff.out',exists_diff

  diff_i = fltarr(n_kinetic,n_field,n_moment,n_r,n_time)
  read_array,diff_i,'diff_i.out',exists_diff_i

  diff_i_trapped = fltarr(n_kinetic,n_field,2,n_r,n_time)
  read_array,diff_i_trapped,'diff_i_trapped.out',exists_diff_i_t

  diff_n = fltarr(n_kinetic,n_field,n_moment,n_n,n_time)
  read_array,diff_n,'diff_n.out',exists_diff_n

  diff_trapped = fltarr(n_kinetic,n_field,2,n_time)
  read_array,diff_trapped,'diff_trapped.out',exists_diff_t

  ;;------------------------------------------------

  ;;--------------------------------------------------------
  ;; Fluxes

  p_moment=4

  gbflux = fltarr(n_kinetic,n_field,p_moment,n_time)
  read_array,gbflux,'gbflux.out',exists_gbflux

  gbflux_trapped = fltarr(n_kinetic,n_field,p_moment,n_time)
  read_array,gbflux_trapped,'gbflux_trapped.out',exists_gbflux_trapped

  gbflux_i = fltarr(n_kinetic,n_field,p_moment,n_r,n_time)
  read_array,gbflux_i,'gbflux_i.out',exists_gbflux_i

  gbflux_i_trapped = fltarr(n_kinetic,n_field,p_moment,n_r,n_time)
  read_array,gbflux_i_trapped,'gbflux_i_trapped.out',exists_gbflux_i_trapped

  gbflux_n = fltarr(n_kinetic,n_field,p_moment,n_n,n_time)
  read_array,gbflux_n,'gbflux_n.out',exists_gbflux_n
  ;;--------------------------------------------------------

  diff_QL_n = fltarr(n_kinetic,n_field,2,n_n,n_time)
  read_array,diff_QL_n,'diff_QL_n.out',exists_diff_QL_n

  phi_squared_QL_n = fltarr(n_n,n_time)
  read_array,phi_squared_QL_n,'phi_squared_QL_n.out',exists_phi_squared_QL_n

  g_squared_QL_n = fltarr(3,n_n,n_time)
  read_array,g_squared_QL_n,'g_squared_QL_n.out',exists_g_squared_QL_n

  Tr_np = fltarr(n_r,n_n,n_time)
  read_array,Tr_np,'nonlinear_transfer_n.out',exists_nonlinear_trasnfer_n

  Eng_np = fltarr(n_r,n_n,n_time)
  read_array,Eng_np,'turbulent_energy_n.out',exists_turbulent_energy_n

  ;;--------------------------------------------
  ;; 0 -> phi_fluxave
  ;; 1 -> a_fluxave
  
  zerobar = fltarr(2,n_r,n_time)
  read_array,zerobar,'zerobar.out',exists_zerobar 
  ;;--------------------------------------------

  source = fltarr(n_kinetic,4,n_r,n_time)
  read_array,source,'source.out',exists_source
  
  ;; Need to DEPRECATE read of zm.out eventually
  zmoment = fltarr(n_r,n_kinetic,2,n_time)
  read_array,zmoment,'zm.out',exists_zmoment
  if (exists_zmoment eq 0) then begin
     zmoment = fltarr(n_r,n_kinetic,2,n_time)
     read_array,zmoment,'moment_zero.out',exists_zmoment 
  endif

  ;; Read various ballooning eigenmode files
  read_ballooning_mode

  ;; Read distribution functions
  read_h

  ;;--------------------------------------------------------
  ;; DATA PREPROCESSING
  ;;--------------------------------------------------------

  midplane_mapper

  no_time:

  ;;--------------------------------------------------------
  ;; Common strings
  ;;
  csa_string     = '!3(c!ds!n/a) t'
  kt_rho_string  = '!3k!d!4h!n q!3!ds!n'
  kr_rho_string  = '!3k!dr!n !4q!3!ds!n'
  chi_string     = '!3(c!ds!n/a)!4q!3!ds!n!u2!n'
  ;;--------------------------------------------------------  

  print,'-------------------------------------'
  print,'|              D O N E              |'
  print,'-------------------------------------'
  
end
