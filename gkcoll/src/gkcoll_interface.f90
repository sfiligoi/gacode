!-------------------------------------------------------------------------
! gkcoll_interface.f90
!
! PURPOSE:
!  Provides interface description for GKCOLL.
!
! CALLING SEQUENCE:
!  call gkcoll_init(...)
!  set gkcoll_*_in variables
!  call gkcoll_run(...)
!  get gkcoll_*_out variables
!-------------------------------------------------------------------------

module gkcoll_interface

  implicit none

  ! Input parameters (set to default values from python/gkcoll_dict.py)
  integer :: gkcoll_n_energy_in = 6
  integer :: gkcoll_n_xi_in     = 17
  integer :: gkcoll_n_theta_in  = 17
  integer :: gkcoll_n_radial_in = 4
  real    :: gkcoll_e_max_in    = 6.0
  real    :: gkcoll_delta_t_in  = 0.01
  real    :: gkcoll_max_time_in = 100.0
  real    :: gkcoll_freq_tol_in = 0.001
  integer :: gkcoll_imp_flag_in = 0
  integer :: gkcoll_restart_write_in = 0
  integer :: gkcoll_restart_mode_in  = 0
  real    :: gkcoll_rupwind_eps_in = 0.0
  integer :: gkcoll_rupwind_n_in = 0
  real    :: gkcoll_tupwind_eps_in = 1.0
  integer :: gkcoll_toroidal_model_in = 1
  integer :: gkcoll_toroidal_num_in  = 1
  real    :: gkcoll_rho_in      = 0.0025
  real    :: gkcoll_k_theta_rho_in  = 0.25
  real    :: gkcoll_rmin_in     = 0.5
  real    :: gkcoll_rmaj_in     = 3.0
  integer :: gkcoll_silent_flag_in = 0
  integer :: gkcoll_equilibrium_model_in = 0
  integer :: gkcoll_collision_model_in   = 4
  integer :: gkcoll_profile_model_in     = 1
  integer :: gkcoll_ipccw_in  = -1
  integer :: gkcoll_btccw_in  = -1
  real    :: gkcoll_te_ade_in = 1.0
  real    :: gkcoll_ne_ade_in = 1.0
  real    :: gkcoll_lambda_debye_in = 0.0
  real    :: gkcoll_profile_lambda_debye_scale_in = 0.0
  integer :: gkcoll_n_species_in = 1
  integer :: gkcoll_z_1_in      = 1
  real    :: gkcoll_mass_1_in   = 1.0
  real    :: gkcoll_dens_1_in   = 1.0
  real    :: gkcoll_temp_1_in   = 1.0
  real    :: gkcoll_dlnndr_1_in = 1.0
  real    :: gkcoll_dlntdr_1_in = 1.0
  real    :: gkcoll_nu_1_in     = 0.1
  integer :: gkcoll_z_2_in      = 1
  real    :: gkcoll_mass_2_in   = 1.0
  real    :: gkcoll_dens_2_in   = 0.0
  real    :: gkcoll_temp_2_in   = 1.0
  real    :: gkcoll_dlnndr_2_in = 1.0
  real    :: gkcoll_dlntdr_2_in = 1.0
  real    :: gkcoll_nu_2_in     = 0.0
  integer :: gkcoll_z_3_in      = 1
  real    :: gkcoll_mass_3_in   = 1.0
  real    :: gkcoll_dens_3_in   = 0.0
  real    :: gkcoll_temp_3_in   = 1.0
  real    :: gkcoll_dlnndr_3_in = 1.0
  real    :: gkcoll_dlntdr_3_in = 1.0
  real    :: gkcoll_nu_3_in     = 0.0
  integer :: gkcoll_z_4_in      = 1
  real    :: gkcoll_mass_4_in   = 1.0
  real    :: gkcoll_dens_4_in   = 0.0
  real    :: gkcoll_temp_4_in   = 1.0
  real    :: gkcoll_dlnndr_4_in = 1.0
  real    :: gkcoll_dlntdr_4_in = 1.0
  real    :: gkcoll_nu_4_in     = 0.0
  integer :: gkcoll_z_5_in      = 1
  real    :: gkcoll_mass_5_in   = 1.0
  real    :: gkcoll_dens_5_in   = 0.0
  real    :: gkcoll_temp_5_in   = 1.0
  real    :: gkcoll_dlnndr_5_in = 1.0
  real    :: gkcoll_dlntdr_5_in = 1.0
  real    :: gkcoll_nu_5_in     = 0.0
  integer :: gkcoll_z_6_in      = 1
  real    :: gkcoll_mass_6_in   = 1.0
  real    :: gkcoll_dens_6_in   = 0.0
  real    :: gkcoll_temp_6_in   = 1.0
  real    :: gkcoll_dlnndr_6_in = 1.0
  real    :: gkcoll_dlntdr_6_in = 1.0
  real    :: gkcoll_nu_6_in     = 0.0
  real    :: gkcoll_q_in        = 2.0
  real    :: gkcoll_shat_in     = 1.0
  real    :: gkcoll_shift_in    = 0.0
  real    :: gkcoll_kappa_in    = 1.0
  real    :: gkcoll_s_kappa_in  = 0.0
  real    :: gkcoll_delta_in    = 0.0
  real    :: gkcoll_s_delta_in  = 0.0
  real    :: gkcoll_zeta_in     = 0.0
  real    :: gkcoll_s_zeta_in   = 0.0
  real    :: gkcoll_zmag_in     = 0.0
  real    :: gkcoll_s_zmag_in   = 0.0
  integer :: gkcoll_geo_ny_in   = 0
  real, dimension(8,0:16) :: gkcoll_geo_yin_in = 0.0

  ! Output parameters
  ! drift-kinetic soln real, dimension(6) :: gkcoll_pflux_dke_out    = 0.0
  real, dimension(6) :: gkcoll_efluxtot_dke_out = 0.0
  real, dimension(6) :: gkcoll_efluxncv_dke_out = 0.0
  real, dimension(6) :: gkcoll_mflux_dke_out    = 0.0
  real, dimension(6) :: gkcoll_vpol_dke_out     = 0.0
  real, dimension(6) :: gkcoll_vtor_dke_out     = 0.0
  real               :: gkcoll_jpar_dke_out     = 0.0
  ! error checking
  integer :: gkcoll_error_status_out=0
  character(len=80) :: gkcoll_error_message_out=''

contains

  ! Map GLOBAL variables to INTERFACE parameters
  subroutine map_global2interface()

    use gkcoll_globals

    implicit none

    gkcoll_n_energy_in = n_energy
    gkcoll_n_xi_in     = n_xi
    gkcoll_n_theta_in  = n_theta
    gkcoll_n_radial_in = n_radial
    gkcoll_e_max_in    = e_max
    gkcoll_delta_t_in  = delta_t
    gkcoll_max_time_in = max_time
    gkcoll_freq_tol_in = freq_tol
    gkcoll_imp_flag_in = imp_flag
    gkcoll_restart_write_in = restart_write
    gkcoll_restart_mode_in  = restart_mode
    gkcoll_rupwind_eps_in = rupwind_eps
    gkcoll_rupwind_n_in = rupwind_n
    gkcoll_tupwind_eps_in = tupwind_eps
    gkcoll_toroidal_model_in = toroidal_model
    gkcoll_toroidal_num_in = toroidal_num
    gkcoll_rho_in      = rho
    gkcoll_k_theta_rho_in  = k_theta_rho
    gkcoll_rmin_in     = rmin
    gkcoll_rmaj_in     = rmaj
    gkcoll_silent_flag_in       = silent_flag
    gkcoll_equilibrium_model_in = equilibrium_model
    gkcoll_collision_model_in   = collision_model
    gkcoll_profile_model_in     = profile_model
    gkcoll_ipccw_in  = ipccw
    gkcoll_btccw_in  = btccw
    gkcoll_te_ade_in = te_ade
    gkcoll_ne_ade_in = ne_ade
    gkcoll_lambda_debye_in = lambda_debye
    gkcoll_profile_lambda_debye_scale_in = profile_lambda_debye_scale
    gkcoll_n_species_in = n_species
    gkcoll_z_1_in      = z(1)
    gkcoll_mass_1_in   = mass(1)
    gkcoll_dens_1_in   = dens(1)
    gkcoll_temp_1_in   = temp(1)
    gkcoll_dlnndr_1_in = dlnndr(1)
    gkcoll_dlntdr_1_in = dlntdr(1)
    gkcoll_nu_1_in     = nu(1)
    gkcoll_z_2_in      = z(2)
    gkcoll_mass_2_in   = mass(2)
    gkcoll_dens_2_in   = dens(2)
    gkcoll_temp_2_in   = temp(2)
    gkcoll_dlnndr_2_in = dlnndr(2)
    gkcoll_dlntdr_2_in = dlntdr(2)
    gkcoll_nu_2_in     = nu(2)
    gkcoll_z_3_in      = z(3)
    gkcoll_mass_3_in   = mass(3)
    gkcoll_dens_3_in   = dens(3)
    gkcoll_temp_3_in   = temp(3)
    gkcoll_dlnndr_3_in = dlnndr(3)
    gkcoll_dlntdr_3_in = dlntdr(3)
    gkcoll_nu_3_in     = nu(3)
    gkcoll_z_4_in      = z(4)
    gkcoll_mass_4_in   = mass(4)
    gkcoll_dens_4_in   = dens(4)
    gkcoll_temp_4_in   = temp(4)
    gkcoll_dlnndr_4_in = dlnndr(4)
    gkcoll_dlntdr_4_in = dlntdr(4)
    gkcoll_nu_4_in     = nu(4)
    gkcoll_z_5_in      = z(5)
    gkcoll_mass_5_in   = mass(5)
    gkcoll_dens_5_in   = dens(5)
    gkcoll_temp_5_in   = temp(5)
    gkcoll_dlnndr_5_in = dlnndr(5)
    gkcoll_dlntdr_5_in = dlntdr(5)
    gkcoll_nu_5_in     = nu(5)
    gkcoll_z_6_in      = z(6)
    gkcoll_mass_6_in   = mass(6)
    gkcoll_dens_6_in   = dens(6)
    gkcoll_temp_6_in   = temp(6)
    gkcoll_dlnndr_6_in = dlnndr(6)
    gkcoll_dlntdr_6_in = dlntdr(6)
    gkcoll_nu_6_in     = nu(6)
    gkcoll_q_in        = q
    gkcoll_shat_in     = shat
    gkcoll_shift_in    = shift
    gkcoll_kappa_in    = kappa
    gkcoll_s_kappa_in  = s_kappa
    gkcoll_delta_in    = delta    
    gkcoll_s_delta_in  = s_delta
    gkcoll_zeta_in     = zeta
    gkcoll_s_zeta_in   = s_zeta
    gkcoll_zmag_in     = zmag
    gkcoll_s_zmag_in   = s_zmag
    gkcoll_geo_ny_in   = geo_ny_in
    gkcoll_geo_yin_in(:,:) = geo_yin_in(:,:)

  end subroutine map_global2interface

  ! Map INTERFACE parameters to GLOBAL variables
  subroutine map_interface2global()

    use gkcoll_globals

    implicit none

    n_energy = gkcoll_n_energy_in
    n_xi     = gkcoll_n_xi_in
    n_theta  = gkcoll_n_theta_in
    n_radial = gkcoll_n_radial_in
    e_max    = gkcoll_e_max_in
    delta_t  = gkcoll_delta_t_in
    max_time = gkcoll_max_time_in
    freq_tol = gkcoll_freq_tol_in
    imp_flag = gkcoll_imp_flag_in
    restart_write = gkcoll_restart_write_in
    restart_mode  = gkcoll_restart_mode_in  
    rupwind_eps = gkcoll_rupwind_eps_in
    rupwind_n = gkcoll_rupwind_n_in
    tupwind_eps = gkcoll_tupwind_eps_in
    toroidal_model = gkcoll_toroidal_model_in
    toroidal_num   = gkcoll_toroidal_num_in
    rho       = gkcoll_rho_in
    k_theta_rho  = gkcoll_k_theta_rho_in
    rmin     = gkcoll_rmin_in
    rmaj     = gkcoll_rmaj_in
    silent_flag       = gkcoll_silent_flag_in
    equilibrium_model = gkcoll_equilibrium_model_in
    collision_model   = gkcoll_collision_model_in
    profile_model     = gkcoll_profile_model_in
    ipccw     = gkcoll_ipccw_in
    btccw     = gkcoll_btccw_in
    te_ade    = gkcoll_te_ade_in
    ne_ade    = gkcoll_ne_ade_in
    lambda_debye = gkcoll_lambda_debye_in
    profile_lambda_debye_scale = gkcoll_profile_lambda_debye_scale_in
    n_species = gkcoll_n_species_in
    z(1)      = gkcoll_z_1_in
    mass(1)   = gkcoll_mass_1_in
    dens(1)   = gkcoll_dens_1_in
    temp(1)   = gkcoll_temp_1_in
    dlnndr(1) = gkcoll_dlnndr_1_in
    dlntdr(1) = gkcoll_dlntdr_1_in
    nu(1)     = gkcoll_nu_1_in
    z(2)      = gkcoll_z_2_in
    mass(2)   = gkcoll_mass_2_in
    dens(2)   = gkcoll_dens_2_in
    temp(2)   = gkcoll_temp_2_in
    dlnndr(2) = gkcoll_dlnndr_2_in
    dlntdr(2) = gkcoll_dlntdr_2_in
    nu(2)     = gkcoll_nu_2_in
    z(3)      = gkcoll_z_3_in
    mass(3)   = gkcoll_mass_3_in
    dens(3)   = gkcoll_dens_3_in
    temp(3)   = gkcoll_temp_3_in
    dlnndr(3) = gkcoll_dlnndr_3_in
    dlntdr(3) = gkcoll_dlntdr_3_in
    nu(3)     = gkcoll_nu_3_in
    z(4)      = gkcoll_z_4_in
    mass(4)   = gkcoll_mass_4_in
    dens(4)   = gkcoll_dens_4_in
    temp(4)   = gkcoll_temp_4_in
    dlnndr(4) = gkcoll_dlnndr_4_in
    dlntdr(4) = gkcoll_dlntdr_4_in
    nu(4)     = gkcoll_nu_4_in
    z(5)      = gkcoll_z_5_in
    mass(5)   = gkcoll_mass_5_in
    dens(5)   = gkcoll_dens_5_in
    temp(5)   = gkcoll_temp_5_in
    dlnndr(5) = gkcoll_dlnndr_5_in
    dlntdr(5) = gkcoll_dlntdr_5_in
    nu(5)     = gkcoll_nu_5_in
    z(6)      = gkcoll_z_6_in
    mass(6)   = gkcoll_mass_6_in
    dens(6)   = gkcoll_dens_6_in
    temp(6)   = gkcoll_temp_6_in
    dlnndr(6) = gkcoll_dlnndr_6_in
    dlntdr(6) = gkcoll_dlntdr_6_in
    nu(6)     = gkcoll_nu_6_in
    q         = gkcoll_q_in
    shat      = gkcoll_shat_in
    shift     = gkcoll_shift_in
    kappa     = gkcoll_kappa_in
    s_kappa   = gkcoll_s_kappa_in
    delta     = gkcoll_delta_in
    s_delta   = gkcoll_s_delta_in
    zeta      = gkcoll_zeta_in
    s_zeta    = gkcoll_s_zeta_in
    zmag      = gkcoll_zmag_in
    s_zmag    = gkcoll_s_zmag_in
    geo_ny_in    = gkcoll_geo_ny_in
    geo_yin_in(:,:) = gkcoll_geo_yin_in(:,:)

  end subroutine map_interface2global

end module gkcoll_interface
