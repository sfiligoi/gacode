!-------------------------------------------------------------------------
! neo_interface.f90
!
! PURPOSE:
!  Provides interface description for NEO.
!
! CALLING SEQUENCE:
!  call neo_init(...)
!  set neo_*_in variables
!  call neo_run(...)
!  get neo_*_out variables
!-------------------------------------------------------------------------

module neo_interface

  implicit none

  ! Input parameters (set to default values from python/neo_dict.py)
  integer :: neo_n_energy_in = 6
  integer :: neo_n_xi_in = 17
  integer :: neo_n_theta_in = 17
  integer :: neo_n_radial_in = 1
  integer :: neo_matsz_scalefac_in = 500
  real    :: neo_rmin_over_a_in = 0.5
  real    :: neo_rmin_over_a_2_in = 0.6
  real    :: neo_rmaj_over_a_in = 3.0
  integer :: neo_silent_flag_in = 0
  integer :: neo_sim_model_in = 1
  integer :: neo_equilibrium_model_in = 0
  integer :: neo_collision_model_in = 4
  integer :: neo_profile_model_in = 1
  integer :: neo_profile_erad0_model_in = 1
  integer :: neo_profile_equilibrium_model_in = 1
  integer :: neo_ipccw_in = -1
  integer :: neo_btccw_in = -1
  real    :: neo_te_ade_in = 1.0
  real    :: neo_ne_ade_in = 1.0
  integer :: neo_rotation_model_in = 1
  real    :: neo_omega_rot_in = 0.0
  real    :: neo_omega_rot_deriv_in = 0.0
  integer :: neo_spitzer_model_in = 0
  real    :: neo_epar0_spitzer_in = 1.0
  integer :: neo_n_species_in = 1
  real    :: neo_nu_1_in = 0.1
  integer :: neo_z_1_in = 1
  real    :: neo_mass_1_in = 1.0
  real    :: neo_dens_1_in = 1.0
  real    :: neo_temp_1_in = 1.0
  real    :: neo_dlnndr_1_in = 1.0
  real    :: neo_dlntdr_1_in = 1.0
  integer :: neo_z_2_in = 1
  real    :: neo_mass_2_in = 1.0
  real    :: neo_dens_2_in = 0.0
  real    :: neo_temp_2_in = 1.0
  real    :: neo_dlnndr_2_in = 1.0
  real    :: neo_dlntdr_2_in = 1.0
  integer :: neo_z_3_in = 1
  real    :: neo_mass_3_in = 1.0
  real    :: neo_dens_3_in = 0.0
  real    :: neo_temp_3_in = 1.0
  real    :: neo_dlnndr_3_in = 1.0
  real    :: neo_dlntdr_3_in = 1.0
  integer :: neo_z_4_in = 1
  real    :: neo_mass_4_in = 1.0
  real    :: neo_dens_4_in = 0.0
  real    :: neo_temp_4_in = 1.0
  real    :: neo_dlnndr_4_in = 1.0
  real    :: neo_dlntdr_4_in = 1.0
  integer :: neo_z_5_in = 1
  real    :: neo_mass_5_in = 1.0
  real    :: neo_dens_5_in = 0.0
  real    :: neo_temp_5_in = 1.0
  real    :: neo_dlnndr_5_in = 1.0
  real    :: neo_dlntdr_5_in = 1.0
  integer :: neo_z_6_in = 1
  real    :: neo_mass_6_in = 1.0
  real    :: neo_dens_6_in = 0.0
  real    :: neo_temp_6_in = 1.0
  real    :: neo_dlnndr_6_in = 1.0
  real    :: neo_dlntdr_6_in = 1.0
  real    :: neo_dphi0dr_in = 0.0
  real    :: neo_epar0_in = 0.0
  real    :: neo_q_in = 2.0
  real    :: neo_rho_star_in = 0.001
  real    :: neo_shear_in = 1.0
  real    :: neo_shift_in = 0.0
  real    :: neo_kappa_in = 1.0
  real    :: neo_s_kappa_in = 0.0
  real    :: neo_delta_in = 0.0
  real    :: neo_s_delta_in = 0.0
  real    :: neo_zeta_in = 0.0
  real    :: neo_s_zeta_in = 0.0
  real    :: neo_zmag_over_a_in = 0.0
  real    :: neo_s_zmag_in = 0.0
  real    :: neo_profile_delta_scale_in = 1.0
  real    :: neo_profile_zeta_scale_in = 1.0
  real    :: neo_profile_zmag_scale_in = 1.0
  integer :: neo_geo_ny_in = 0
  real, dimension(8,0:16) :: neo_geo_yin_in = 0.0

  ! Output parameters
  ! theory
  real    :: neo_pflux_thHH_out  = 0.0
  real    :: neo_eflux_thHHi_out = 0.0
  real    :: neo_eflux_thHHe_out = 0.0
  real    :: neo_eflux_thCHi_out = 0.0
  real, dimension(6) :: neo_pflux_thHS_out = 0.0
  real, dimension(6) :: neo_eflux_thHS_out = 0.0
  real    :: neo_jpar_thHH_out   = 0.0
  real    :: neo_jpar_thS_out    = 0.0
  ! drift-kinetic soln
  real, dimension(6) :: neo_pflux_dke_out    = 0.0
  real, dimension(6) :: neo_efluxtot_dke_out = 0.0
  real, dimension(6) :: neo_efluxncv_dke_out = 0.0
  real, dimension(6) :: neo_mflux_dke_out    = 0.0
  real, dimension(6) :: neo_vpol_dke_out     = 0.0
  real, dimension(6) :: neo_vtor_dke_out     = 0.0
  real               :: neo_jpar_dke_out     = 0.0
  ! gyro-viscosity
  real, dimension(6) :: neo_pflux_gv_out  = 0.0
  real, dimension(6) :: neo_efluxtot_gv_out  = 0.0
  real, dimension(6) :: neo_efluxncv_gv_out  = 0.0
  real, dimension(6) :: neo_mflux_gv_out  = 0.0
  ! error checking
  integer :: neo_error_status_out=0
  character(len=80) :: neo_error_message_out=''

contains

  ! Map GLOBAL variables to INTERFACE parameters
  subroutine map_global2interface()

    use neo_globals

    implicit none

    neo_n_energy_in = n_energy
    neo_n_xi_in = n_xi
    neo_n_theta_in = n_theta
    neo_n_radial_in = n_radial
    neo_matsz_scalefac_in = matsz_scalefac
    neo_rmin_over_a_in = rmin_1_in
    neo_rmin_over_a_in = rmin_1_in
    neo_rmin_over_a_2_in = rmin_2_in
    neo_rmaj_over_a_in = rmaj_in
    neo_silent_flag_in = silent_flag
    neo_sim_model_in = sim_model
    neo_equilibrium_model_in = equilibrium_model
    neo_collision_model_in = collision_model
    neo_profile_model_in = profile_model
    neo_profile_erad0_model_in = profile_erad0_model
    neo_profile_equilibrium_model_in = profile_equilibrium_model
    neo_ipccw_in  = ipccw_in
    neo_btccw_in  = btccw_in
    neo_te_ade_in = te_ade_in
    neo_ne_ade_in = ne_ade_in
    neo_rotation_model_in = rotation_model
    neo_omega_rot_in = omega_rot_in
    neo_omega_rot_deriv_in = omega_rot_deriv_in
    neo_spitzer_model_in = spitzer_model
    neo_epar0_spitzer_in = epar0_spitzer
    neo_n_species_in = n_species
    neo_nu_1_in = nu_1_in
    neo_z_1_in = z_in(1)
    neo_mass_1_in = mass_in(1)
    neo_dens_1_in = dens_in(1)
    neo_temp_1_in = temp_in(1)
    neo_dlnndr_1_in = dlnndr_in(1)
    neo_dlntdr_1_in = dlntdr_in(1)
    neo_z_2_in = z_in(2)
    neo_mass_2_in = mass_in(2)
    neo_dens_2_in = dens_in(2)
    neo_temp_2_in = temp_in(2)
    neo_dlnndr_2_in = dlnndr_in(2)
    neo_dlntdr_2_in = dlntdr_in(2)
    neo_z_3_in = z_in(3)
    neo_mass_3_in = mass_in(3)
    neo_dens_3_in = dens_in(3)
    neo_temp_3_in = temp_in(3)
    neo_dlnndr_3_in = dlnndr_in(3)
    neo_dlntdr_3_in = dlntdr_in(3)
    neo_z_4_in = z_in(4)
    neo_mass_4_in = mass_in(4)
    neo_dens_4_in = dens_in(4)
    neo_temp_4_in = temp_in(4)
    neo_dlnndr_4_in = dlnndr_in(4)
    neo_dlntdr_4_in = dlntdr_in(4)
    neo_z_5_in = z_in(5)
    neo_mass_5_in = mass_in(5)
    neo_dens_5_in = dens_in(5)
    neo_temp_5_in = temp_in(5)
    neo_dlnndr_5_in = dlnndr_in(5)
    neo_dlntdr_5_in = dlntdr_in(5)
    neo_z_6_in = z_in(6)
    neo_mass_6_in = mass_in(6)
    neo_dens_6_in = dens_in(6)
    neo_temp_6_in = temp_in(6)
    neo_dlnndr_6_in = dlnndr_in(6)
    neo_dlntdr_6_in = dlntdr_in(6)
    neo_dphi0dr_in = dphi0dr_in
    neo_epar0_in = epar0_in
    neo_q_in = q_in
    neo_rho_star_in = rho_in
    neo_shear_in = shat_in
    neo_shift_in = shift_in
    neo_kappa_in = kappa_in
    neo_s_kappa_in = s_kappa_in
    neo_delta_in = delta_in
    neo_s_delta_in = s_delta_in
    neo_zeta_in = zeta_in
    neo_s_zeta_in = s_zeta_in
    neo_zmag_over_a_in = zmag_in
    neo_s_zmag_in = s_zmag_in
    neo_profile_delta_scale_in = profile_delta_scale
    neo_profile_zeta_scale_in = profile_zeta_scale
    neo_profile_zmag_scale_in = profile_zmag_scale
    neo_geo_ny_in = geo_ny_in
    neo_geo_yin_in(:,:) = geo_yin_in(:,:)

  end subroutine map_global2interface

  ! Map INTERFACE parameters to GLOBAL variables
  subroutine map_interface2global()

    use neo_globals

    implicit none

    n_energy = neo_n_energy_in
    n_xi = neo_n_xi_in
    n_theta = neo_n_theta_in
    n_radial = neo_n_radial_in
    matsz_scalefac = neo_matsz_scalefac_in
    rmin_1_in = neo_rmin_over_a_in
    rmin_2_in = neo_rmin_over_a_2_in
    rmaj_in = neo_rmaj_over_a_in
    silent_flag = neo_silent_flag_in
    sim_model = neo_sim_model_in
    equilibrium_model = neo_equilibrium_model_in
    collision_model = neo_collision_model_in
    profile_model = neo_profile_model_in
    profile_erad0_model = neo_profile_erad0_model_in
    profile_equilibrium_model = neo_profile_equilibrium_model_in
    ipccw_in  = neo_ipccw_in
    btccw_in  = neo_btccw_in
    te_ade_in = neo_te_ade_in
    ne_ade_in = neo_ne_ade_in
    rotation_model = neo_rotation_model_in
    omega_rot_in = neo_omega_rot_in
    omega_rot_deriv_in = neo_omega_rot_deriv_in
    spitzer_model = neo_spitzer_model_in
    epar0_spitzer = neo_epar0_spitzer_in
    n_species = neo_n_species_in
    nu_1_in = neo_nu_1_in
    z_in(1) = neo_z_1_in
    mass_in(1) = neo_mass_1_in
    dens_in(1) = neo_dens_1_in
    temp_in(1) = neo_temp_1_in
    dlnndr_in(1) = neo_dlnndr_1_in
    dlntdr_in(1) = neo_dlntdr_1_in
    z_in(2) = neo_z_2_in
    mass_in(2) = neo_mass_2_in
    dens_in(2) = neo_dens_2_in
    temp_in(2) = neo_temp_2_in
    dlnndr_in(2) = neo_dlnndr_2_in
    dlntdr_in(2) = neo_dlntdr_2_in
    z_in(3) = neo_z_3_in
    mass_in(3) = neo_mass_3_in
    dens_in(3) = neo_dens_3_in
    temp_in(3) = neo_temp_3_in
    dlnndr_in(3) = neo_dlnndr_3_in
    dlntdr_in(3) = neo_dlntdr_3_in
    z_in(4) = neo_z_4_in
    mass_in(4) = neo_mass_4_in
    dens_in(4) = neo_dens_4_in
    temp_in(4) = neo_temp_4_in
    dlnndr_in(4) = neo_dlnndr_4_in
    dlntdr_in(4) = neo_dlntdr_4_in
    z_in(5) = neo_z_5_in
    mass_in(5) = neo_mass_5_in
    dens_in(5) = neo_dens_5_in
    temp_in(5) = neo_temp_5_in
    dlnndr_in(5) = neo_dlnndr_5_in
    dlntdr_in(5) = neo_dlntdr_5_in
    z_in(6) = neo_z_6_in
    mass_in(6) = neo_mass_6_in
    dens_in(6) = neo_dens_6_in
    temp_in(6) = neo_temp_6_in
    dlnndr_in(6) = neo_dlnndr_6_in
    dlntdr_in(6) = neo_dlntdr_6_in
    dphi0dr_in = neo_dphi0dr_in
    epar0_in = neo_epar0_in
    q_in = neo_q_in
    rho_in = neo_rho_star_in
    shat_in = neo_shear_in
    shift_in = neo_shift_in
    kappa_in = neo_kappa_in
    s_kappa_in = neo_s_kappa_in
    delta_in = neo_delta_in
    s_delta_in = neo_s_delta_in
    zeta_in = neo_zeta_in
    s_zeta_in = neo_s_zeta_in
    zmag_in = neo_zmag_over_a_in
    s_zmag_in = neo_s_zmag_in
    profile_delta_scale = neo_profile_delta_scale_in
    profile_zeta_scale = neo_profile_zeta_scale_in
    profile_zmag_scale = neo_profile_zmag_scale_in
    geo_ny_in = neo_geo_ny_in
    geo_yin_in(:,:) = neo_geo_yin_in(:,:)

  end subroutine map_interface2global

end module neo_interface
