!-------------------------------------------------------------------------
! cgyro_interface.f90
!
! PURPOSE:
!  Provides interface description for CGYRO
!
! CALLING SEQUENCE:
!  call cgyro_init(...)
!  set cgyro_*_in variables
!  call cgyro_run(...)
!  get cgyro_*_out variables
!-------------------------------------------------------------------------

module cgyro_interface

  implicit none

  ! Input parameters (set to default values from python/cgyro_dict.py)
  integer :: cgyro_n_energy_in = 6
  integer :: cgyro_n_xi_in     = 16
  integer :: cgyro_n_theta_in  = 31
  integer :: cgyro_n_radial_in = 4
  integer :: cgyro_e_max_in    = 7
  real    :: cgyro_delta_t_in  = 0.01
  real    :: cgyro_max_time_in = 100.0
  real    :: cgyro_freq_tol_in = 0.001
  integer :: cgyro_restart_write_in = 0
  integer :: cgyro_restart_mode_in  = 0
  real    :: cgyro_up_radial_in = 1.0
  integer :: cgyro_up_radial_n_in = 4
  real    :: cgyro_tupwind_eps_in = 1.0
  integer :: cgyro_toroidal_model_in = 1
  integer :: cgyro_toroidal_num_in  = 1
  real    :: cgyro_rho_in      = 0.0025
  real    :: cgyro_k_theta_rho_in  = 0.25
  real    :: cgyro_r_length_rho_in  = 64.0
  real    :: cgyro_rmin_in     = 0.5
  real    :: cgyro_rmaj_in     = 3.0
  integer :: cgyro_silent_flag_in = 0
  integer :: cgyro_equilibrium_model_in = 2
  integer :: cgyro_collision_model_in   = 1
  real    :: cgyro_te_ade_in = 1.0
  real    :: cgyro_ne_ade_in = 1.0
  real    :: cgyro_lambda_debye_in = 0.0
  integer :: cgyro_n_species_in = 1
  integer :: cgyro_z_1_in      = 1
  real    :: cgyro_mass_1_in   = 1.0
  real    :: cgyro_dens_1_in   = 1.0
  real    :: cgyro_temp_1_in   = 1.0
  real    :: cgyro_dlnndr_1_in = 1.0
  real    :: cgyro_dlntdr_1_in = 1.0
  real    :: cgyro_nu_1_in     = 0.1
  integer :: cgyro_z_2_in      = 1
  real    :: cgyro_mass_2_in   = 1.0
  real    :: cgyro_dens_2_in   = 0.0
  real    :: cgyro_temp_2_in   = 1.0
  real    :: cgyro_dlnndr_2_in = 1.0
  real    :: cgyro_dlntdr_2_in = 1.0
  real    :: cgyro_nu_2_in     = 0.0
  integer :: cgyro_z_3_in      = 1
  real    :: cgyro_mass_3_in   = 1.0
  real    :: cgyro_dens_3_in   = 0.0
  real    :: cgyro_temp_3_in   = 1.0
  real    :: cgyro_dlnndr_3_in = 1.0
  real    :: cgyro_dlntdr_3_in = 1.0
  real    :: cgyro_nu_3_in     = 0.0
  integer :: cgyro_z_4_in      = 1
  real    :: cgyro_mass_4_in   = 1.0
  real    :: cgyro_dens_4_in   = 0.0
  real    :: cgyro_temp_4_in   = 1.0
  real    :: cgyro_dlnndr_4_in = 1.0
  real    :: cgyro_dlntdr_4_in = 1.0
  real    :: cgyro_nu_4_in     = 0.0
  integer :: cgyro_z_5_in      = 1
  real    :: cgyro_mass_5_in   = 1.0
  real    :: cgyro_dens_5_in   = 0.0
  real    :: cgyro_temp_5_in   = 1.0
  real    :: cgyro_dlnndr_5_in = 1.0
  real    :: cgyro_dlntdr_5_in = 1.0
  real    :: cgyro_nu_5_in     = 0.0
  integer :: cgyro_z_6_in      = 1
  real    :: cgyro_mass_6_in   = 1.0
  real    :: cgyro_dens_6_in   = 0.0
  real    :: cgyro_temp_6_in   = 1.0
  real    :: cgyro_dlnndr_6_in = 1.0
  real    :: cgyro_dlntdr_6_in = 1.0
  real    :: cgyro_nu_6_in     = 0.0
  real    :: cgyro_q_in        = 2.0
  real    :: cgyro_s_in     = 1.0
  real    :: cgyro_shift_in    = 0.0
  real    :: cgyro_kappa_in    = 1.0
  real    :: cgyro_s_kappa_in  = 0.0
  real    :: cgyro_delta_in    = 0.0
  real    :: cgyro_s_delta_in  = 0.0
  real    :: cgyro_zeta_in     = 0.0
  real    :: cgyro_s_zeta_in   = 0.0
  real    :: cgyro_zmag_in     = 0.0
  real    :: cgyro_s_zmag_in   = 0.0
  real    :: cgyro_beta_star_in  = 0.0
  integer :: cgyro_n_field_in    = 1
  real    :: cgyro_betae_unit_in = 0.0
  integer :: cgyro_geo_ny_in   = 0
  real, dimension(8,0:16) :: cgyro_geo_yin_in = 0.0

  ! Output parameters
  ! drift-kinetic soln real, dimension(6) :: cgyro_pflux_dke_out    = 0.0
  real, dimension(6) :: cgyro_efluxtot_dke_out = 0.0
  real, dimension(6) :: cgyro_efluxncv_dke_out = 0.0
  real, dimension(6) :: cgyro_mflux_dke_out    = 0.0
  real, dimension(6) :: cgyro_vpol_dke_out     = 0.0
  real, dimension(6) :: cgyro_vtor_dke_out     = 0.0
  real               :: cgyro_jpar_dke_out     = 0.0
  ! error checking
  integer :: cgyro_error_status_out=0
  character(len=80) :: cgyro_error_message_out=''

contains

  ! Map GLOBAL variables to INTERFACE parameters
  subroutine map_global2interface()

    use cgyro_globals

    implicit none

    cgyro_n_energy_in = n_energy
    cgyro_n_xi_in     = n_xi
    cgyro_n_theta_in  = n_theta
    cgyro_n_radial_in = n_radial
    cgyro_e_max_in    = e_max
    cgyro_delta_t_in  = delta_t
    cgyro_max_time_in = max_time
    cgyro_freq_tol_in = freq_tol
    cgyro_restart_write_in = restart_write
    cgyro_restart_mode_in  = restart_mode
    cgyro_up_radial_in = up_radial
    cgyro_up_radial_n_in = up_radial_n
    cgyro_tupwind_eps_in = tupwind_eps
    cgyro_toroidal_model_in = toroidal_model
    cgyro_toroidal_num_in = toroidal_num
    cgyro_rho_in      = rho
    cgyro_k_theta_rho_in  = k_theta_rho
    cgyro_r_length_rho_in = r_length_rho
    cgyro_rmin_in     = rmin
    cgyro_rmaj_in     = rmaj
    cgyro_silent_flag_in       = silent_flag
    cgyro_equilibrium_model_in = equilibrium_model
    cgyro_collision_model_in   = collision_model
    cgyro_te_ade_in = te_ade
    cgyro_ne_ade_in = ne_ade
    cgyro_lambda_debye_in = lambda_debye
    cgyro_n_species_in = n_species
    cgyro_z_1_in      = z(1)
    cgyro_mass_1_in   = mass(1)
    cgyro_dens_1_in   = dens(1)
    cgyro_temp_1_in   = temp(1)
    cgyro_dlnndr_1_in = dlnndr(1)
    cgyro_dlntdr_1_in = dlntdr(1)
    cgyro_nu_1_in     = nu(1)
    cgyro_z_2_in      = z(2)
    cgyro_mass_2_in   = mass(2)
    cgyro_dens_2_in   = dens(2)
    cgyro_temp_2_in   = temp(2)
    cgyro_dlnndr_2_in = dlnndr(2)
    cgyro_dlntdr_2_in = dlntdr(2)
    cgyro_nu_2_in     = nu(2)
    cgyro_z_3_in      = z(3)
    cgyro_mass_3_in   = mass(3)
    cgyro_dens_3_in   = dens(3)
    cgyro_temp_3_in   = temp(3)
    cgyro_dlnndr_3_in = dlnndr(3)
    cgyro_dlntdr_3_in = dlntdr(3)
    cgyro_nu_3_in     = nu(3)
    cgyro_z_4_in      = z(4)
    cgyro_mass_4_in   = mass(4)
    cgyro_dens_4_in   = dens(4)
    cgyro_temp_4_in   = temp(4)
    cgyro_dlnndr_4_in = dlnndr(4)
    cgyro_dlntdr_4_in = dlntdr(4)
    cgyro_nu_4_in     = nu(4)
    cgyro_z_5_in      = z(5)
    cgyro_mass_5_in   = mass(5)
    cgyro_dens_5_in   = dens(5)
    cgyro_temp_5_in   = temp(5)
    cgyro_dlnndr_5_in = dlnndr(5)
    cgyro_dlntdr_5_in = dlntdr(5)
    cgyro_nu_5_in     = nu(5)
    cgyro_z_6_in      = z(6)
    cgyro_mass_6_in   = mass(6)
    cgyro_dens_6_in   = dens(6)
    cgyro_temp_6_in   = temp(6)
    cgyro_dlnndr_6_in = dlnndr(6)
    cgyro_dlntdr_6_in = dlntdr(6)
    cgyro_nu_6_in     = nu(6)
    cgyro_q_in        = q
    cgyro_s_in     = s
    cgyro_shift_in    = shift
    cgyro_kappa_in    = kappa
    cgyro_s_kappa_in  = s_kappa
    cgyro_delta_in    = delta    
    cgyro_s_delta_in  = s_delta
    cgyro_zeta_in     = zeta
    cgyro_s_zeta_in   = s_zeta
    cgyro_zmag_in     = zmag
    cgyro_s_zmag_in   = s_zmag
    cgyro_beta_star_in  = beta_star
    cgyro_n_field_in    = n_field
    cgyro_betae_unit_in = betae_unit
    cgyro_geo_ny_in   = geo_ny_in
    cgyro_geo_yin_in(:,:) = geo_yin_in(:,:)

  end subroutine map_global2interface

  ! Map INTERFACE parameters to GLOBAL variables
  subroutine map_interface2global()

    use cgyro_globals

    implicit none

    n_energy = cgyro_n_energy_in
    n_xi     = cgyro_n_xi_in
    n_theta  = cgyro_n_theta_in
    n_radial = cgyro_n_radial_in
    e_max    = cgyro_e_max_in
    delta_t  = cgyro_delta_t_in
    max_time = cgyro_max_time_in
    freq_tol = cgyro_freq_tol_in
    restart_write = cgyro_restart_write_in
    restart_mode  = cgyro_restart_mode_in  
    up_radial = cgyro_up_radial_in
    up_radial_n = cgyro_up_radial_n_in
    tupwind_eps = cgyro_tupwind_eps_in
    toroidal_model = cgyro_toroidal_model_in
    toroidal_num   = cgyro_toroidal_num_in
    rho       = cgyro_rho_in
    k_theta_rho  = cgyro_k_theta_rho_in
    r_length_rho = cgyro_r_length_rho_in
    rmin     = cgyro_rmin_in
    rmaj     = cgyro_rmaj_in
    silent_flag       = cgyro_silent_flag_in
    equilibrium_model = cgyro_equilibrium_model_in
    collision_model   = cgyro_collision_model_in
    te_ade    = cgyro_te_ade_in
    ne_ade    = cgyro_ne_ade_in
    lambda_debye = cgyro_lambda_debye_in
    n_species = cgyro_n_species_in
    z(1)      = cgyro_z_1_in
    mass(1)   = cgyro_mass_1_in
    dens(1)   = cgyro_dens_1_in
    temp(1)   = cgyro_temp_1_in
    dlnndr(1) = cgyro_dlnndr_1_in
    dlntdr(1) = cgyro_dlntdr_1_in
    nu(1)     = cgyro_nu_1_in
    z(2)      = cgyro_z_2_in
    mass(2)   = cgyro_mass_2_in
    dens(2)   = cgyro_dens_2_in
    temp(2)   = cgyro_temp_2_in
    dlnndr(2) = cgyro_dlnndr_2_in
    dlntdr(2) = cgyro_dlntdr_2_in
    nu(2)     = cgyro_nu_2_in
    z(3)      = cgyro_z_3_in
    mass(3)   = cgyro_mass_3_in
    dens(3)   = cgyro_dens_3_in
    temp(3)   = cgyro_temp_3_in
    dlnndr(3) = cgyro_dlnndr_3_in
    dlntdr(3) = cgyro_dlntdr_3_in
    nu(3)     = cgyro_nu_3_in
    z(4)      = cgyro_z_4_in
    mass(4)   = cgyro_mass_4_in
    dens(4)   = cgyro_dens_4_in
    temp(4)   = cgyro_temp_4_in
    dlnndr(4) = cgyro_dlnndr_4_in
    dlntdr(4) = cgyro_dlntdr_4_in
    nu(4)     = cgyro_nu_4_in
    z(5)      = cgyro_z_5_in
    mass(5)   = cgyro_mass_5_in
    dens(5)   = cgyro_dens_5_in
    temp(5)   = cgyro_temp_5_in
    dlnndr(5) = cgyro_dlnndr_5_in
    dlntdr(5) = cgyro_dlntdr_5_in
    nu(5)     = cgyro_nu_5_in
    z(6)      = cgyro_z_6_in
    mass(6)   = cgyro_mass_6_in
    dens(6)   = cgyro_dens_6_in
    temp(6)   = cgyro_temp_6_in
    dlnndr(6) = cgyro_dlnndr_6_in
    dlntdr(6) = cgyro_dlntdr_6_in
    nu(6)     = cgyro_nu_6_in
    q         = cgyro_q_in
    s      = cgyro_s_in
    shift     = cgyro_shift_in
    kappa     = cgyro_kappa_in
    s_kappa   = cgyro_s_kappa_in
    delta     = cgyro_delta_in
    s_delta   = cgyro_s_delta_in
    zeta      = cgyro_zeta_in
    s_zeta    = cgyro_s_zeta_in
    zmag      = cgyro_zmag_in
    s_zmag    = cgyro_s_zmag_in
    beta_star  = cgyro_beta_star_in
    n_field    = cgyro_n_field_in
    betae_unit = cgyro_betae_unit_in
    geo_ny_in  = cgyro_geo_ny_in
    geo_yin_in(:,:) = cgyro_geo_yin_in(:,:)

  end subroutine map_interface2global

end module cgyro_interface
