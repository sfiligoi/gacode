subroutine tgyro_neo_map

  use tgyro_globals
  use neo_interface

  implicit none
  real :: mu1
  real :: gamma_p0,u000
  integer :: is

  ! Initialize NEO
  call neo_init(paths(i_r-1),gyro_comm)

  neo_test_flag_in = gyrotest_flag
  
  ! Simulation mode (dke solve vs. analytic)
  if (loc_neo_method == 1) then
     ! Analytic theory
     neo_sim_model_in = 0
  else
     ! Kinetic NEO and theory (no NCLASS)
     neo_sim_model_in = 2
  end if

  ! Resolution 
  neo_n_energy_in = 5
  neo_n_xi_in = 11
  neo_n_theta_in = tgyro_neo_n_theta

  ! Geometry
  neo_equilibrium_model_in = 2
  neo_rmin_over_a_in       = r(i_r)/r_min
  neo_rmaj_over_a_in       = r_maj(i_r)/r_min
  neo_q_in                 = abs(q(i_r))
  neo_shear_in             = s(i_r)
  neo_shift_in             = shift(i_r)
  neo_kappa_in             = kappa(i_r)
  neo_s_kappa_in           = s_kappa(i_r)
  neo_delta_in             = delta(i_r)
  neo_s_delta_in           = s_delta(i_r)
  neo_zeta_in              = zeta(i_r)
  neo_s_zeta_in            = s_zeta(i_r)
  neo_zmag_over_a_in       = zmag(i_r)
  neo_s_zmag_in            = dzmag(i_r)

  neo_ipccw_in = -signb*signq
  neo_btccw_in = -signb

  neo_rho_star_in  = 0.001

  neo_n_species_in = loc_n_ion+1
  if (loc_n_ion > 9) then
     call tgyro_catch_error('ERROR: (TGYRO) n_ion > 9 not supported in NEO interface.') 
  endif

  ! Assuming NEO m_norm = mi(1) (first ion)
  mu1 = sqrt(mi(1)/(me*loc_me_multiplier))
  
  ! Electrons
  neo_z_in(2)      = -1
  neo_mass_in(2)   = 1.0/mu1**2
  neo_dens_in(2)   = 1.0
  neo_temp_in(2)   = 1.0/(ti(1,i_r)/te(i_r))
  neo_dlnndr_in(2) = r_min*dlnnedr(i_r)
  neo_dlntdr_in(2) = r_min*dlntedr(i_r)

  ! Main ions
  neo_z_in(1)      = int(zi_vec(1))
  neo_mass_in(1)   = 1.0
  neo_dens_in(1)   = ni(1,i_r)/ne(i_r)
  neo_temp_in(1)   = 1.0
  neo_dlnndr_in(1) = r_min*dlnnidr(1,i_r)
  neo_dlntdr_in(1) = r_min*dlntidr(1,i_r)
  neo_nu_1_in      = nui(1,i_r)*r_min/c_s(i_r)/sqrt(ti(1,i_r)/te(i_r))

  ! Ion #2
  do is=2,neo_n_species_in

     neo_z_in(is+1)      = int(zi_vec(is))
     neo_mass_in(is+1)   = 1.0/sqrt(mi(1)/mi(is))**2
     neo_dens_in(is+1)   = ni(is,i_r)/ne(i_r)
     neo_temp_in(is+1)   = ti(is,i_r)/ti(1,i_r)
     neo_dlnndr_in(is+1) = r_min*dlnnidr(is,i_r)
     neo_dlntdr_in(is+1) = r_min*dlntidr(is,i_r)

  enddo

  ! Rotation is always *on* in NEO.
  !
  ! COORDINATES: The signs of all rotation-related quantities below are 
  ! inherited (unchanged) from input.profiles.  In general NEO expects 
  ! these to be correctly signed/oriented.

  gamma_p0  = -r_maj(i_r)*w0p(i_r)
  u000      = r_maj(i_r)*w0(i_r)

  neo_rotation_model_in = 2
  neo_omega_rot_in = u000 * r_min / r_maj(i_r) &
       / (c_s(i_r) * sqrt(ti(1,i_r)/te(i_r)))
  neo_omega_rot_deriv_in = -gamma_p0 * r_min**2 / r_maj(i_r) &
       / (c_s(i_r) * sqrt(ti(1,i_r)/te(i_r)))

  ! Parameter only used for global runs.
  neo_rmin_over_a_2_in = neo_rmin_over_a_in

  ! General geometry Fourier coefficients

  if (loc_num_equil_flag == 1) then
     ! Resetting if numerical equilibrium wanted.
     neo_equilibrium_model_in = 3
  endif

  neo_geo_ny_in = n_fourier_geo
  neo_geo_yin_in(:,:) = a_fourier_geo(:,:,i_r)

end subroutine tgyro_neo_map
