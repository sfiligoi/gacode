!------------------------------------------------------------
! tgyro_glf23_map.f90
!
! PURPOSE:
!  Mapping from TGYRO internal variables to GLF23 interface.
!------------------------------------------------------------

subroutine tgyro_glf23_map

  use tgyro_globals
  use glf23_interface

  implicit none

  ! Local variables
  integer :: i_ion
  real :: q_abs
  real :: gamma_eb0
  real :: gamma_p0

  q_abs = abs(q(i_r))

  ! Initialize GLF23
  call glf23_init(paths(i_r-1),gyro_comm)

  ! put_model_parameters 
  !----------------------------------------------------------------
  ! Want fluxes from GLF23
  glf23_use_transport_model_in = .true.
  !
  glf23_version_in = tgyro_glf23_revision
  !
  !----------------------------------------------------------------
  ! put_species 
  !
  ! Number of species (max=3)
  !
  glf23_ns_in = loc_n_ion+1
  if (glf23_ns_in > 3) then
     call tgyro_catch_error('ERROR: (tgyro_glf23_map) Too many ions in GLF23.')
  endif
  ! Charges: e,i,z
  glf23_zs_in(1) = -1.0
  do i_ion=1,loc_n_ion
     glf23_zs_in(i_ion+1) = zi_vec(i_ion)
  enddo
  ! Assuming mi(1) is normalizing mass
  ! Mass ratios: me/mi(1),mi(1)/mi(1),m(2)/mi(1),...
  glf23_mass_in(1) = (me*loc_me_multiplier)/mi(1)
  do i_ion=1,loc_n_ion 
     glf23_mass_in(i_ion+1) = mi(i_ion)/mi(1)
  enddo
  !
  !----------------------------------------------------------------
  ! put_gradients
  !-----------------------------------
  ! Density gradients (e,i,z)
  glf23_rlns_in(1) = r_min*dlnnedr(i_r)
  do i_ion=1,loc_n_ion
     glf23_rlns_in(i_ion+1) = r_min*dlnnidr(i_ion,i_r)
  enddo
  ! Temperature gradients (e,i,z)
  glf23_rlts_in(1) = r_min*dlntedr(i_r)
  do i_ion=1,loc_n_ion
     glf23_rlts_in(i_ion+1) = r_min*dlntidr(i_ion,i_r)
  enddo
  ! Gamma_ExB (ExB shearing rate, units of a/cs)
  if (tgyro_rotation_flag == 1) then
     gamma_p0  = -r_maj(i_r)*f_rot(i_r)*w0_norm
     gamma_eb0 = gamma_p0*r(i_r)/(q_abs*r_maj(i_r))
  ! sign of vexb_shear does not matter for GLF23 
     glf23_vexb_shear_in    = gamma_eb0*r_min/c_s(i_r)  
     glf23_vpar_shear_in(:) = gamma_p0*r_min/c_s(i_r)
  endif
  !-----------------------------------
  ! put_averages
  !----------------------------------------------------------------
  ! Te/Te,Ti/Te,Tz/Te
  glf23_taus_in(1) = 1.0
  do i_ion=1,loc_n_ion
     glf23_taus_in(i_ion+1) = ti(i_ion,i_r)/te(i_r)
  enddo
  ! ne/ne,ni(1)/ne,ni(2)/ne,...
  glf23_as_in(1) = 1.0
  do i_ion=1,loc_n_ion
     glf23_as_in(i_ion+1) = ni(i_ion,i_r)/ne(i_r) 
  enddo
  ! Electron beta used for electromagnetic calculations
  glf23_betae_in = betae_unit(i_r)*loc_betae_scale
  ! Electron collision frequency
  glf23_xnue_in = nue(i_r)*r_min/c_s(i_r)*loc_nu_scale
  !
  !----------------------------------------------------------------
  !put_s_alpha_geometry
  !
  glf23_rmin_sa_in     = r(i_r)/r_min
  glf23_rmaj_sa_in     = r_maj(i_r)/r_min
  glf23_q_sa_in        = q_abs
  glf23_shat_sa_in     = s(i_r)
  glf23_alpha_sa_in    = r_maj(i_r)*beta_unit(i_r)*dlnpdr(i_r)*q_abs**2
  glf23_xwell_sa_in    = 0.0
  glf23_theta0_sa_in   = 0.0
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! glf23_run output control 
  if (tgyro_glf23_dump_flag == 0) then
     glf23_dump_flag_in   = .false.
  else
     glf23_dump_flag_in = .true.
  endif
!
end subroutine tgyro_glf23_map
