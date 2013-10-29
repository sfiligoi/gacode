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
  integer :: neo_matsz_scalefac_in = 50
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
  real    :: neo_dlntdre_ade_in = 1.0
  real    :: neo_dlnndre_ade_in = 1.0
  integer :: neo_rotation_model_in = 1
  real    :: neo_omega_rot_in = 0.0
  real    :: neo_omega_rot_deriv_in = 0.0
  integer :: neo_spitzer_model_in = 0
  real    :: neo_epar0_spitzer_in = 1.0
  integer :: neo_coll_uncoupledei_model_in = 0
  integer :: neo_n_species_in = 1
  real    :: neo_nu_1_in = 0.1
  integer :: neo_z_1_in = 1
  real    :: neo_mass_1_in = 1.0
  real    :: neo_dens_1_in = 1.0
  real    :: neo_temp_1_in = 1.0
  real    :: neo_dlnndr_1_in = 1.0
  real    :: neo_dlntdr_1_in = 1.0
  real    :: neo_profile_dlnndr_1_scale_in = 1.0
  real    :: neo_profile_dlntdr_1_scale_in = 1.0
  integer :: neo_z_2_in = 1
  real    :: neo_mass_2_in = 1.0
  real    :: neo_dens_2_in = 0.0
  real    :: neo_temp_2_in = 1.0
  real    :: neo_dlnndr_2_in = 1.0
  real    :: neo_dlntdr_2_in = 1.0
  real    :: neo_profile_dlnndr_2_scale_in = 1.0
  real    :: neo_profile_dlntdr_2_scale_in = 1.0
  integer :: neo_z_3_in = 1
  real    :: neo_mass_3_in = 1.0
  real    :: neo_dens_3_in = 0.0
  real    :: neo_temp_3_in = 1.0
  real    :: neo_dlnndr_3_in = 1.0
  real    :: neo_dlntdr_3_in = 1.0
  real    :: neo_profile_dlnndr_3_scale_in = 1.0
  real    :: neo_profile_dlntdr_3_scale_in = 1.0
  integer :: neo_z_4_in = 1
  real    :: neo_mass_4_in = 1.0
  real    :: neo_dens_4_in = 0.0
  real    :: neo_temp_4_in = 1.0
  real    :: neo_dlnndr_4_in = 1.0
  real    :: neo_dlntdr_4_in = 1.0
  real    :: neo_profile_dlnndr_4_scale_in = 1.0
  real    :: neo_profile_dlntdr_4_scale_in = 1.0
  integer :: neo_z_5_in = 1
  real    :: neo_mass_5_in = 1.0
  real    :: neo_dens_5_in = 0.0
  real    :: neo_temp_5_in = 1.0
  real    :: neo_dlnndr_5_in = 1.0
  real    :: neo_dlntdr_5_in = 1.0
  real    :: neo_profile_dlnndr_5_scale_in = 1.0
  real    :: neo_profile_dlntdr_5_scale_in = 1.0
  integer :: neo_z_6_in = 1
  real    :: neo_mass_6_in = 1.0
  real    :: neo_dens_6_in = 0.0
  real    :: neo_temp_6_in = 1.0
  real    :: neo_dlnndr_6_in = 1.0
  real    :: neo_dlntdr_6_in = 1.0
  real    :: neo_profile_dlnndr_6_scale_in = 1.0
  real    :: neo_profile_dlntdr_6_scale_in = 1.0
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
  integer :: neo_threed_model_in = 0
  integer :: neo_threed_exb_model_in = 0
  real    :: neo_threed_exb_dphi0dr_in = 0.0
  integer :: neo_geo_ny_in = 0
  real, dimension(8,0:32) :: neo_geo_yin_in = 0.0
  ! the exception of the default is subroutine_flag
  integer :: neo_subroutine_flag = 1

  ! Output parameters
  ! theory
  real    :: neo_pflux_thHH_out  = 0.0
  real    :: neo_eflux_thHHi_out = 0.0
  real    :: neo_eflux_thHHe_out = 0.0
  real    :: neo_eflux_thCHi_out = 0.0
  real, dimension(6) :: neo_pflux_thHS_out = 0.0
  real, dimension(6) :: neo_eflux_thHS_out = 0.0
  real    :: neo_jpar_thS_out    = 0.0
  real    :: neo_jpar_thK_out    = 0.0
  real    :: neo_jpar_thN_out    = 0.0
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
  ! nclass viscosities
  real, dimension(6) :: neo_nclassvis_out  = 0.0
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
    neo_dlntdre_ade_in = dlntdre_ade_in
    neo_dlnndre_ade_in = dlnndre_ade_in
    neo_rotation_model_in = rotation_model
    neo_omega_rot_in = omega_rot_in
    neo_omega_rot_deriv_in = omega_rot_deriv_in
    neo_spitzer_model_in = spitzer_model
    neo_epar0_spitzer_in = epar0_spitzer
    neo_coll_uncoupledei_model_in = coll_uncoupledei_model
    neo_n_species_in = n_species
    neo_nu_1_in = nu_1_in
    neo_z_1_in = z_in(1)
    neo_mass_1_in = mass_in(1)
    neo_dens_1_in = dens_in(1)
    neo_temp_1_in = temp_in(1)
    neo_dlnndr_1_in = dlnndr_in(1)
    neo_dlntdr_1_in = dlntdr_in(1)
    neo_profile_dlnndr_1_scale_in = profile_dlnndr_scale(1)
    neo_profile_dlntdr_1_scale_in = profile_dlntdr_scale(1)
    neo_z_2_in = z_in(2)
    neo_mass_2_in = mass_in(2)
    neo_dens_2_in = dens_in(2)
    neo_temp_2_in = temp_in(2)
    neo_dlnndr_2_in = dlnndr_in(2)
    neo_dlntdr_2_in = dlntdr_in(2)
    neo_profile_dlnndr_2_scale_in = profile_dlnndr_scale(2)
    neo_profile_dlntdr_2_scale_in = profile_dlntdr_scale(2)
    neo_z_3_in = z_in(3)
    neo_mass_3_in = mass_in(3)
    neo_dens_3_in = dens_in(3)
    neo_temp_3_in = temp_in(3)
    neo_dlnndr_3_in = dlnndr_in(3)
    neo_dlntdr_3_in = dlntdr_in(3)
    neo_profile_dlnndr_3_scale_in = profile_dlnndr_scale(3)
    neo_profile_dlntdr_3_scale_in = profile_dlntdr_scale(3)
    neo_z_4_in = z_in(4)
    neo_mass_4_in = mass_in(4)
    neo_dens_4_in = dens_in(4)
    neo_temp_4_in = temp_in(4)
    neo_dlnndr_4_in = dlnndr_in(4)
    neo_dlntdr_4_in = dlntdr_in(4)
    neo_profile_dlnndr_4_scale_in = profile_dlnndr_scale(4)
    neo_profile_dlntdr_4_scale_in = profile_dlntdr_scale(4)
    neo_z_5_in = z_in(5)
    neo_mass_5_in = mass_in(5)
    neo_dens_5_in = dens_in(5)
    neo_temp_5_in = temp_in(5)
    neo_dlnndr_5_in = dlnndr_in(5)
    neo_dlntdr_5_in = dlntdr_in(5)
    neo_profile_dlnndr_5_scale_in = profile_dlnndr_scale(5)
    neo_profile_dlntdr_5_scale_in = profile_dlntdr_scale(5)
    neo_z_6_in = z_in(6)
    neo_mass_6_in = mass_in(6)
    neo_dens_6_in = dens_in(6)
    neo_temp_6_in = temp_in(6)
    neo_dlnndr_6_in = dlnndr_in(6)
    neo_dlntdr_6_in = dlntdr_in(6)
    neo_profile_dlnndr_6_scale_in = profile_dlnndr_scale(6)
    neo_profile_dlntdr_6_scale_in = profile_dlntdr_scale(6)
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
    neo_threed_model_in = threed_model
    neo_threed_exb_model_in = threed_exb_model
    neo_threed_exb_dphi0dr_in = threed_exb_dphi0dr
    neo_geo_ny_in = geo_ny_in
    neo_geo_yin_in(:,:) = geo_yin_in(:,:)
    neo_subroutine_flag = subroutine_flag

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
    dlntdre_ade_in = neo_dlntdre_ade_in
    dlnndre_ade_in = neo_dlnndre_ade_in
    rotation_model = neo_rotation_model_in
    omega_rot_in = neo_omega_rot_in
    omega_rot_deriv_in = neo_omega_rot_deriv_in
    spitzer_model = neo_spitzer_model_in
    epar0_spitzer = neo_epar0_spitzer_in
    coll_uncoupledei_model = neo_coll_uncoupledei_model_in
    n_species = neo_n_species_in
    nu_1_in = neo_nu_1_in
    z_in(1) = neo_z_1_in
    mass_in(1) = neo_mass_1_in
    dens_in(1) = neo_dens_1_in
    temp_in(1) = neo_temp_1_in
    dlnndr_in(1) = neo_dlnndr_1_in
    dlntdr_in(1) = neo_dlntdr_1_in
    profile_dlnndr_scale(1) = neo_profile_dlnndr_1_scale_in
    profile_dlntdr_scale(1) = neo_profile_dlntdr_1_scale_in
    z_in(2) = neo_z_2_in
    mass_in(2) = neo_mass_2_in
    dens_in(2) = neo_dens_2_in
    temp_in(2) = neo_temp_2_in
    dlnndr_in(2) = neo_dlnndr_2_in
    dlntdr_in(2) = neo_dlntdr_2_in
    profile_dlnndr_scale(2) = neo_profile_dlnndr_2_scale_in
    profile_dlntdr_scale(2) = neo_profile_dlntdr_2_scale_in
    z_in(3) = neo_z_3_in
    mass_in(3) = neo_mass_3_in
    dens_in(3) = neo_dens_3_in
    temp_in(3) = neo_temp_3_in
    dlnndr_in(3) = neo_dlnndr_3_in
    dlntdr_in(3) = neo_dlntdr_3_in
    profile_dlnndr_scale(3) = neo_profile_dlnndr_3_scale_in
    profile_dlntdr_scale(3) = neo_profile_dlntdr_3_scale_in
    z_in(4) = neo_z_4_in
    mass_in(4) = neo_mass_4_in
    dens_in(4) = neo_dens_4_in
    temp_in(4) = neo_temp_4_in
    dlnndr_in(4) = neo_dlnndr_4_in
    dlntdr_in(4) = neo_dlntdr_4_in
    profile_dlnndr_scale(4) = neo_profile_dlnndr_4_scale_in
    profile_dlntdr_scale(4) = neo_profile_dlntdr_4_scale_in
    z_in(5) = neo_z_5_in
    mass_in(5) = neo_mass_5_in
    dens_in(5) = neo_dens_5_in
    temp_in(5) = neo_temp_5_in
    dlnndr_in(5) = neo_dlnndr_5_in
    dlntdr_in(5) = neo_dlntdr_5_in
    profile_dlnndr_scale(5) = neo_profile_dlnndr_5_scale_in
    profile_dlntdr_scale(5) = neo_profile_dlntdr_5_scale_in
    z_in(6) = neo_z_6_in
    mass_in(6) = neo_mass_6_in
    dens_in(6) = neo_dens_6_in
    temp_in(6) = neo_temp_6_in
    dlnndr_in(6) = neo_dlnndr_6_in
    dlntdr_in(6) = neo_dlntdr_6_in
    profile_dlnndr_scale(6) = neo_profile_dlnndr_6_scale_in
    profile_dlntdr_scale(6) = neo_profile_dlntdr_6_scale_in
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
    threed_model = neo_threed_model_in
    threed_exb_model = neo_threed_exb_model_in
    threed_exb_dphi0dr = neo_threed_exb_dphi0dr_in
    geo_ny_in = neo_geo_ny_in
    geo_yin_in(:,:) = neo_geo_yin_in(:,:)
    subroutine_flag = neo_subroutine_flag

    call interfacelocaldump

  end subroutine map_interface2global

  subroutine interfacelocaldump
    use neo_globals
    implicit none
    
    if(neo_silent_flag_in > 0 .or. i_proc > 0) return

    open(unit=1,file=trim(path)//'out.neo.localdump',status='replace')
    
    if(neo_n_radial_in > 1) then
       write(1,*) 'Localdump only works with n_radial=1'
       close(1)
       return
    endif
    if(neo_profile_model_in /= 1) then
       write(1,*) 'Localdump only works with profile_model=1'
       close(1)
       return
    endif
    if(neo_equilibrium_model_in == 3) then
       write(1,*) 'Localdump does not work with general geometry'
       close(1)
       return
    endif
    
    write(1,20) 'N_ENERGY=',neo_n_energy_in
    write(1,20) 'N_XI=',neo_n_xi_in
    write(1,20) 'N_THETA=',neo_n_theta_in
    write(1,20) 'N_RADIAL=',neo_n_radial_in
    write(1,20) 'MATSZ_SCALEFAC=',neo_matsz_scalefac_in
    write(1,30) 'RMIN_OVER_A=',neo_rmin_over_a_in
    write(1,30) 'RMAJ_OVER_A=',neo_rmaj_over_a_in
    write(1,20) 'SILENT_FLAG=',neo_silent_flag_in
    write(1,20) 'SIM_MODEL=',neo_sim_model_in
    write(1,20) 'EQUILIBRIUM_MODEL=',neo_equilibrium_model_in
    write(1,20) 'COLLISION_MODEL=',neo_collision_model_in
    write(1,20) 'IPCCW=',neo_ipccw_in
    write(1,20) 'BTCCW=',neo_btccw_in
    write(1,30) 'TE_ADE=',neo_te_ade_in
    write(1,30) 'NE_ADE=',neo_ne_ade_in
    write(1,30) 'DLNTDRE_ADE=',neo_dlntdre_ade_in
    write(1,30) 'DLNNDRE_ADE=',neo_dlnndre_ade_in
    write(1,20) 'ROTATION_MODEL=',neo_rotation_model_in
    write(1,30) 'OMEGA_ROT=',neo_omega_rot_in
    write(1,30) 'OMEGA_ROT_DERIV=',neo_omega_rot_deriv_in
    write(1,20) 'SPITZER_MODEL=',neo_spitzer_model_in
    write(1,30) 'EPAR0_SPITZER=',neo_epar0_spitzer_in
    write(1,20) 'COLL_UNCOUPLEDEI_MODEL=',neo_coll_uncoupledei_model_in
    write(1,20) 'N_SPECIES=',neo_n_species_in
    write(1,30) 'NU_1=',neo_nu_1_in
    
    write(1,20) 'Z_1=',neo_z_1_in
    write(1,30) 'MASS_1=',neo_mass_1_in
    write(1,30) 'DENS_1=',neo_dens_1_in
    write(1,30) 'TEMP_1=',neo_temp_1_in
    write(1,30) 'DLNNDR_1=',neo_dlnndr_1_in
    write(1,30) 'DLNTDR_1=',neo_dlntdr_1_in
    if (neo_n_species_in > 1) then
       write(1,20) 'Z_2=',neo_z_2_in
       write(1,30) 'MASS_2=',neo_mass_2_in
       write(1,30) 'DENS_2=',neo_dens_2_in
       write(1,30) 'TEMP_2=',neo_temp_2_in
       write(1,30) 'DLNNDR_2=',neo_dlnndr_2_in
       write(1,30) 'DLNTDR_2=',neo_dlntdr_2_in
    endif
    if (neo_n_species_in > 2) then
       write(1,20) 'Z_3=',neo_z_3_in
       write(1,30) 'MASS_3=',neo_mass_3_in
       write(1,30) 'DENS_3=',neo_dens_3_in
       write(1,30) 'TEMP_3=',neo_temp_3_in
       write(1,30) 'DLNNDR_3=',neo_dlnndr_3_in
       write(1,30) 'DLNTDR_3=',neo_dlntdr_3_in
    endif
    if (neo_n_species_in > 3) then
       write(1,20) 'Z_4=',neo_z_4_in
       write(1,30) 'MASS_4=',neo_mass_4_in
       write(1,30) 'DENS_4=',neo_dens_4_in
       write(1,30) 'TEMP_4=',neo_temp_4_in
       write(1,30) 'DLNNDR_4=',neo_dlnndr_4_in
       write(1,30) 'DLNTDR_4=',neo_dlntdr_4_in
    endif
    if (neo_n_species_in > 4) then
       write(1,20) 'Z_5=',neo_z_5_in
       write(1,30) 'MASS_5=',neo_mass_5_in
       write(1,30) 'DENS_5=',neo_dens_5_in
       write(1,30) 'TEMP_5=',neo_temp_5_in
       write(1,30) 'DLNNDR_5=',neo_dlnndr_5_in
       write(1,30) 'DLNTDR_5=',neo_dlntdr_5_in
    endif
    if (neo_n_species_in > 5) then
       write(1,20) 'Z_6=',neo_z_6_in
       write(1,30) 'MASS_6=',neo_mass_6_in
       write(1,30) 'DENS_6=',neo_dens_6_in
       write(1,30) 'TEMP_6=',neo_temp_6_in
       write(1,30) 'DLNNDR_6=',neo_dlnndr_6_in
       write(1,30) 'DLNTDR_6=',neo_dlntdr_6_in
    endif
  
    write(1,30) 'DPHI0DR=',neo_dphi0dr_in
    write(1,30) 'EPAR0=',neo_epar0_in
    write(1,30) 'Q=',neo_q_in
    write(1,30) 'RHO_STAR=',neo_rho_star_in
    write(1,30) 'SHEAR=',neo_shear_in
    write(1,30) 'SHIFT=',neo_shift_in
    write(1,30) 'KAPPA=',neo_kappa_in
    write(1,30) 'S_KAPPA=',neo_s_kappa_in
    write(1,30) 'DELTA=',neo_delta_in
    write(1,30) 'S_DELTA=',neo_s_delta_in
    write(1,30) 'ZETA=',neo_zeta_in
    write(1,30) 'S_ZETA=',neo_s_zeta_in
    write(1,30) 'ZMAG_OVER_A=',neo_zmag_over_a_in
    write(1,30) 'S_ZMAG=',neo_s_zmag_in
    
    close(1)
    
20  format(t2,a,i3)
30  format(t2,a,1pe12.5)
    
  end subroutine interfacelocaldump
  
  
end module neo_interface
