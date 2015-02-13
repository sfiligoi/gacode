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

  ! Input parameters 
  ! 
  ! Defaults should agree with values in cgyro_parse.py
  !
  integer :: cgyro_n_energy_in = 6
  integer :: cgyro_n_xi_in     = 16
  integer :: cgyro_n_theta_in  = 16
  integer :: cgyro_n_radial_in = 4
  integer :: cgyro_n_toroidal_in = 1
  integer :: cgyro_n_field_in = 1
  integer :: cgyro_e_max_in    = 7
  real    :: cgyro_delta_t_in  = 0.01
  real    :: cgyro_max_time_in = 10.0
  integer :: cgyro_print_step_in = 50
  real    :: cgyro_freq_tol_in = 0.005
  integer :: cgyro_restart_write_in = 1
  integer :: cgyro_restart_mode_in  = 0
  real    :: cgyro_up_radial_in = 1.0
  integer :: cgyro_up_radial_n_in = 4
  real    :: cgyro_up_theta_in = 1.0
  real    :: cgyro_up_xi_in = 1.0
  integer :: cgyro_implicit_flag_in = 0
  real    :: cgyro_ky_in = 0.3
  real    :: cgyro_box_size_in  = 1.0
  integer :: cgyro_silent_flag_in = 0
  integer :: cgyro_equilibrium_model_in = 2
  integer :: cgyro_collision_model_in   = 4
  integer :: cgyro_collision_mom_restore_in = 1
  integer :: cgyro_collision_ene_restore_in = 1
  integer :: cgyro_collision_ene_diffusion_in = 1
  integer :: cgyro_collision_kperp_in = 1
  integer :: cgyro_collision_field_model_in = 1
  integer :: cgyro_collision_trap_model_in = 1
  integer :: cgyro_zf_test_flag_in = 0
  integer :: cgyro_nonlinear_flag_in = 0
  real    :: cgyro_te_ade_in = 1.0
  real    :: cgyro_ne_ade_in = 1.0
  real    :: cgyro_masse_ade_in = 0.0002724486
  real    :: cgyro_lambda_debye_in = 0.0

  real    :: cgyro_rmin_in    = 0.5
  real    :: cgyro_rmaj_in    = 3.0
  real    :: cgyro_q_in       = 2.0
  real    :: cgyro_s_in       = 1.0
  real    :: cgyro_shift_in   = 0.0
  real    :: cgyro_kappa_in   = 1.0
  real    :: cgyro_s_kappa_in = 0.0
  real    :: cgyro_delta_in   = 0.0
  real    :: cgyro_s_delta_in = 0.0
  real    :: cgyro_zeta_in    = 0.0
  real    :: cgyro_s_zeta_in  = 0.0
  real    :: cgyro_zmag_in    = 0.0
  real    :: cgyro_s_zmag_in  = 0.0
  real    :: cgyro_beta_star_in  = 0.0
  real    :: cgyro_betae_unit_in = 0.0

  integer :: cgyro_n_species_in = 1
  real    :: cgyro_nu_ee_in      = 0.1
  integer, dimension(6) :: cgyro_z_in = 1
  real, dimension(6)    :: cgyro_mass_in = 1.0
  real, dimension(6)    :: cgyro_dens_in = 1.0
  real, dimension(6)    :: cgyro_temp_in = 1.0
  real, dimension(6)    :: cgyro_dlnndr_in = 1.0
  real, dimension(6)    :: cgyro_dlntdr_in = 1.0

  ! Fourier geometry
  integer :: cgyro_geo_ny_in   = 0
  real, dimension(8,0:16) :: cgyro_geo_yin_in = 0.0

  ! Error checking
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
    cgyro_n_toroidal_in = n_toroidal
    cgyro_e_max_in    = e_max
    cgyro_delta_t_in  = delta_t
    cgyro_max_time_in = max_time
    cgyro_freq_tol_in = freq_tol
    cgyro_restart_write_in = restart_write
    cgyro_restart_mode_in  = restart_mode
    cgyro_up_radial_in = up_radial
    cgyro_up_radial_n_in = up_radial_n
    cgyro_up_theta_in = up_theta
    cgyro_up_xi_in = up_xi
    cgyro_implicit_flag_in = implicit_flag
    cgyro_rmin_in     = rmin
    cgyro_rmaj_in     = rmaj
    cgyro_silent_flag_in       = silent_flag
    cgyro_equilibrium_model_in = equilibrium_model
    cgyro_collision_model_in   = collision_model
    cgyro_collision_mom_restore_in = collision_mom_restore
    cgyro_collision_ene_restore_in = collision_ene_restore
    cgyro_collision_ene_diffusion_in = collision_ene_diffusion
    cgyro_collision_kperp_in = collision_kperp
    cgyro_collision_field_model_in = collision_field_model
    cgyro_collision_trap_model_in = collision_trap_model
    cgyro_te_ade_in = te_ade
    cgyro_ne_ade_in = ne_ade
    cgyro_masse_ade_in = masse_ade
    cgyro_lambda_debye_in = lambda_debye
    cgyro_n_species_in = n_species
    cgyro_z_in(:)    = z(:)
    cgyro_mass_in(:) = mass(:)
    cgyro_dens_in(:)   = dens(:)
    cgyro_temp_in(:)   = temp(:)
    cgyro_dlnndr_in(:) = dlnndr(:)
    cgyro_dlntdr_in(:) = dlntdr(:)
    cgyro_q_in         = q
    cgyro_s_in        = s
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
    up_theta = cgyro_up_theta_in
    up_xi = cgyro_up_xi_in
    implicit_flag = cgyro_implicit_flag_in
    rmin     = cgyro_rmin_in
    rmaj     = cgyro_rmaj_in
    silent_flag       = cgyro_silent_flag_in
    equilibrium_model = cgyro_equilibrium_model_in
    collision_model   = cgyro_collision_model_in
    collision_mom_restore = cgyro_collision_mom_restore_in
    collision_ene_restore = cgyro_collision_ene_restore_in
    collision_ene_diffusion = cgyro_collision_ene_diffusion_in
    collision_kperp = cgyro_collision_kperp_in
    collision_field_model = cgyro_collision_field_model_in
    collision_trap_model = cgyro_collision_trap_model_in
    te_ade    = cgyro_te_ade_in
    ne_ade    = cgyro_ne_ade_in
    masse_ade = cgyro_masse_ade_in
    lambda_debye = cgyro_lambda_debye_in
    n_species = cgyro_n_species_in
    z(:)      = cgyro_z_in(:)
    mass(:)   = cgyro_mass_in(:)
    dens(:)   = cgyro_dens_in(:)
    temp(:)   = cgyro_temp_in(:)
    dlnndr(:) = cgyro_dlnndr_in(:)
    dlntdr(:) = cgyro_dlntdr_in(:)
    q         = cgyro_q_in
    s         = cgyro_s_in
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
