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
  integer :: gkcoll_n_xi_in = 17
  integer :: gkcoll_n_theta_in = 17
  integer :: gkcoll_n_radial_in = 1
  real    :: gkcoll_e_max_in    = 6.0
  integer :: gkcoll_matsz_scalefac_in = 500
  real    :: gkcoll_rmin_over_a_in = 0.5
  real    :: gkcoll_rmaj_over_a_in = 3.0
  integer :: gkcoll_silent_flag_in = 0
  integer :: gkcoll_equilibrium_model_in = 0
  integer :: gkcoll_collision_model_in = 4
  integer :: gkcoll_profile_model_in = 1
  integer :: gkcoll_profile_temprescale_model_in = 0
  integer :: gkcoll_profile_equilibrium_model_in = 1
  integer :: gkcoll_ipccw_in = -1
  integer :: gkcoll_btccw_in = -1
  real    :: gkcoll_te_ade_in = 1.0
  real    :: gkcoll_ne_ade_in = 1.0
  integer :: gkcoll_n_species_in = 1
  integer :: gkcoll_z_1_in = 1
  real    :: gkcoll_mass_1_in = 1.0
  real    :: gkcoll_dens_1_in = 1.0
  real    :: gkcoll_temp_1_in = 1.0
  real    :: gkcoll_dlnndr_1_in = 1.0
  real    :: gkcoll_dlntdr_1_in = 1.0
  real    :: gkcoll_nu_1_in = 0.1
  integer :: gkcoll_z_2_in = 1
  real    :: gkcoll_mass_2_in = 1.0
  real    :: gkcoll_dens_2_in = 0.0
  real    :: gkcoll_temp_2_in = 1.0
  real    :: gkcoll_dlnndr_2_in = 1.0
  real    :: gkcoll_dlntdr_2_in = 1.0
  real    :: gkcoll_nu_2_in = 0.0
  integer :: gkcoll_z_3_in = 1
  real    :: gkcoll_mass_3_in = 1.0
  real    :: gkcoll_dens_3_in = 0.0
  real    :: gkcoll_temp_3_in = 1.0
  real    :: gkcoll_dlnndr_3_in = 1.0
  real    :: gkcoll_dlntdr_3_in = 1.0
  real    :: gkcoll_nu_3_in = 0.0
  integer :: gkcoll_z_4_in = 1
  real    :: gkcoll_mass_4_in = 1.0
  real    :: gkcoll_dens_4_in = 0.0
  real    :: gkcoll_temp_4_in = 1.0
  real    :: gkcoll_dlnndr_4_in = 1.0
  real    :: gkcoll_dlntdr_4_in = 1.0
  real    :: gkcoll_nu_4_in = 0.0
  integer :: gkcoll_z_5_in = 1
  real    :: gkcoll_mass_5_in = 1.0
  real    :: gkcoll_dens_5_in = 0.0
  real    :: gkcoll_temp_5_in = 1.0
  real    :: gkcoll_dlnndr_5_in = 1.0
  real    :: gkcoll_dlntdr_5_in = 1.0
  real    :: gkcoll_nu_5_in = 0.0
  integer :: gkcoll_z_6_in = 1
  real    :: gkcoll_mass_6_in = 1.0
  real    :: gkcoll_dens_6_in = 0.0
  real    :: gkcoll_temp_6_in = 1.0
  real    :: gkcoll_dlnndr_6_in = 1.0
  real    :: gkcoll_dlntdr_6_in = 1.0
  real    :: gkcoll_nu_6_in = 0.0
  real    :: gkcoll_q_in = 2.0
  real    :: gkcoll_rho_star_in = 0.001
  real    :: gkcoll_shear_in = 1.0
  real    :: gkcoll_shift_in = 0.0
  real    :: gkcoll_kappa_in = 1.0
  real    :: gkcoll_s_kappa_in = 0.0
  real    :: gkcoll_delta_in = 0.0
  real    :: gkcoll_s_delta_in = 0.0
  real    :: gkcoll_zeta_in = 0.0
  real    :: gkcoll_s_zeta_in = 0.0
  real    :: gkcoll_zmag_over_a_in = 0.0
  real    :: gkcoll_s_zmag_in = 0.0
  real    :: gkcoll_profile_delta_scale_in = 1.0
  real    :: gkcoll_profile_zeta_scale_in = 1.0
  real    :: gkcoll_profile_zmag_scale_in = 1.0
  integer :: gkcoll_geo_ny_in = 0
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
    gkcoll_n_xi_in = n_xi
    gkcoll_n_theta_in = n_theta
    gkcoll_n_radial_in = n_radial
    gkcoll_e_max_in = e_max
    gkcoll_matsz_scalefac_in = matsz_scalefac
    gkcoll_rmin_over_a_in = rmin_in
    gkcoll_rmaj_over_a_in = rmaj_in
    gkcoll_silent_flag_in = silent_flag
    gkcoll_equilibrium_model_in = equilibrium_model
    gkcoll_collision_model_in = collision_model
    gkcoll_profile_model_in = profile_model
    gkcoll_profile_temprescale_model_in = profile_temprescale_model
    gkcoll_profile_equilibrium_model_in = profile_equilibrium_model
    gkcoll_ipccw_in  = ipccw_in
    gkcoll_btccw_in  = btccw_in
    gkcoll_te_ade_in = te_ade_in
    gkcoll_ne_ade_in = ne_ade_in
    gkcoll_n_species_in = n_species
    gkcoll_z_1_in = z_in(1)
    gkcoll_mass_1_in = mass_in(1)
    gkcoll_dens_1_in = dens_in(1)
    gkcoll_temp_1_in = temp_in(1)
    gkcoll_dlnndr_1_in = dlnndr_in(1)
    gkcoll_dlntdr_1_in = dlntdr_in(1)
    gkcoll_nu_1_in = nu_in(1)
    gkcoll_z_2_in = z_in(2)
    gkcoll_mass_2_in = mass_in(2)
    gkcoll_dens_2_in = dens_in(2)
    gkcoll_temp_2_in = temp_in(2)
    gkcoll_dlnndr_2_in = dlnndr_in(2)
    gkcoll_dlntdr_2_in = dlntdr_in(2)
    gkcoll_nu_2_in = nu_in(2)
    gkcoll_z_3_in = z_in(3)
    gkcoll_mass_3_in = mass_in(3)
    gkcoll_dens_3_in = dens_in(3)
    gkcoll_temp_3_in = temp_in(3)
    gkcoll_dlnndr_3_in = dlnndr_in(3)
    gkcoll_dlntdr_3_in = dlntdr_in(3)
    gkcoll_nu_3_in = nu_in(3)
    gkcoll_z_4_in = z_in(4)
    gkcoll_mass_4_in = mass_in(4)
    gkcoll_dens_4_in = dens_in(4)
    gkcoll_temp_4_in = temp_in(4)
    gkcoll_dlnndr_4_in = dlnndr_in(4)
    gkcoll_dlntdr_4_in = dlntdr_in(4)
    gkcoll_nu_4_in = nu_in(4)
    gkcoll_z_5_in = z_in(5)
    gkcoll_mass_5_in = mass_in(5)
    gkcoll_dens_5_in = dens_in(5)
    gkcoll_temp_5_in = temp_in(5)
    gkcoll_dlnndr_5_in = dlnndr_in(5)
    gkcoll_dlntdr_5_in = dlntdr_in(5)
    gkcoll_nu_5_in = nu_in(5)
    gkcoll_z_6_in = z_in(6)
    gkcoll_mass_6_in = mass_in(6)
    gkcoll_dens_6_in = dens_in(6)
    gkcoll_temp_6_in = temp_in(6)
    gkcoll_dlnndr_6_in = dlnndr_in(6)
    gkcoll_dlntdr_6_in = dlntdr_in(6)
    gkcoll_nu_6_in = nu_in(6)
    gkcoll_q_in = q_in
    gkcoll_rho_star_in = rho_in
    gkcoll_shear_in = shat_in
    gkcoll_shift_in = shift_in
    gkcoll_kappa_in = kappa_in
    gkcoll_s_kappa_in = s_kappa_in
    gkcoll_delta_in = delta_in
    gkcoll_s_delta_in = s_delta_in
    gkcoll_zeta_in = zeta_in
    gkcoll_s_zeta_in = s_zeta_in
    gkcoll_zmag_over_a_in = zmag_in
    gkcoll_s_zmag_in = s_zmag_in
    gkcoll_profile_delta_scale_in = profile_delta_scale
    gkcoll_profile_zeta_scale_in = profile_zeta_scale
    gkcoll_profile_zmag_scale_in = profile_zmag_scale
    gkcoll_geo_ny_in = geo_ny_in
    gkcoll_geo_yin_in(:,:) = geo_yin_in(:,:)

  end subroutine map_global2interface

  ! Map INTERFACE parameters to GLOBAL variables
  subroutine map_interface2global()

    use gkcoll_globals

    implicit none

    n_energy = gkcoll_n_energy_in
    n_xi = gkcoll_n_xi_in
    n_theta = gkcoll_n_theta_in
    n_radial = gkcoll_n_radial_in
    e_max = gkcoll_e_max_in
    matsz_scalefac = gkcoll_matsz_scalefac_in
    rmin_in = gkcoll_rmin_over_a_in
    rmaj_in = gkcoll_rmaj_over_a_in
    silent_flag = gkcoll_silent_flag_in
    equilibrium_model = gkcoll_equilibrium_model_in
    collision_model = gkcoll_collision_model_in
    profile_model = gkcoll_profile_model_in
    profile_temprescale_model = gkcoll_profile_temprescale_model_in
    profile_equilibrium_model = gkcoll_profile_equilibrium_model_in
    ipccw_in  = gkcoll_ipccw_in
    btccw_in  = gkcoll_btccw_in
    te_ade_in = gkcoll_te_ade_in
    ne_ade_in = gkcoll_ne_ade_in
    n_species = gkcoll_n_species_in
    z_in(1) = gkcoll_z_1_in
    mass_in(1) = gkcoll_mass_1_in
    dens_in(1) = gkcoll_dens_1_in
    temp_in(1) = gkcoll_temp_1_in
    dlnndr_in(1) = gkcoll_dlnndr_1_in
    dlntdr_in(1) = gkcoll_dlntdr_1_in
    nu_in(1) = gkcoll_nu_1_in
    z_in(2) = gkcoll_z_2_in
    mass_in(2) = gkcoll_mass_2_in
    dens_in(2) = gkcoll_dens_2_in
    temp_in(2) = gkcoll_temp_2_in
    dlnndr_in(2) = gkcoll_dlnndr_2_in
    dlntdr_in(2) = gkcoll_dlntdr_2_in
    nu_in(2) = gkcoll_nu_2_in
    z_in(3) = gkcoll_z_3_in
    mass_in(3) = gkcoll_mass_3_in
    dens_in(3) = gkcoll_dens_3_in
    temp_in(3) = gkcoll_temp_3_in
    dlnndr_in(3) = gkcoll_dlnndr_3_in
    dlntdr_in(3) = gkcoll_dlntdr_3_in
    nu_in(3) = gkcoll_nu_3_in
    z_in(4) = gkcoll_z_4_in
    mass_in(4) = gkcoll_mass_4_in
    dens_in(4) = gkcoll_dens_4_in
    temp_in(4) = gkcoll_temp_4_in
    dlnndr_in(4) = gkcoll_dlnndr_4_in
    dlntdr_in(4) = gkcoll_dlntdr_4_in
    nu_in(4) = gkcoll_nu_4_in
    z_in(5) = gkcoll_z_5_in
    mass_in(5) = gkcoll_mass_5_in
    dens_in(5) = gkcoll_dens_5_in
    temp_in(5) = gkcoll_temp_5_in
    dlnndr_in(5) = gkcoll_dlnndr_5_in
    dlntdr_in(5) = gkcoll_dlntdr_5_in
    nu_in(5) = gkcoll_nu_5_in
    z_in(6) = gkcoll_z_6_in
    mass_in(6) = gkcoll_mass_6_in
    dens_in(6) = gkcoll_dens_6_in
    temp_in(6) = gkcoll_temp_6_in
    dlnndr_in(6) = gkcoll_dlnndr_6_in
    dlntdr_in(6) = gkcoll_dlntdr_6_in
    nu_in(6) = gkcoll_nu_6_in
    q_in = gkcoll_q_in
    rho_in = gkcoll_rho_star_in
    shat_in = gkcoll_shear_in
    shift_in = gkcoll_shift_in
    kappa_in = gkcoll_kappa_in
    s_kappa_in = gkcoll_s_kappa_in
    delta_in = gkcoll_delta_in
    s_delta_in = gkcoll_s_delta_in
    zeta_in = gkcoll_zeta_in
    s_zeta_in = gkcoll_s_zeta_in
    zmag_in = gkcoll_zmag_over_a_in
    s_zmag_in = gkcoll_s_zmag_in
    profile_delta_scale = gkcoll_profile_delta_scale_in
    profile_zeta_scale = gkcoll_profile_zeta_scale_in
    profile_zmag_scale = gkcoll_profile_zmag_scale_in
    geo_ny_in = gkcoll_geo_ny_in
    geo_yin_in(:,:) = gkcoll_geo_yin_in(:,:)

  end subroutine map_interface2global

end module gkcoll_interface
