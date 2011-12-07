subroutine tgyro_neo_map

  use tgyro_globals
  use neo_interface

  implicit none
  real :: mu1

  if (loc_n_ion > 5) then
     call tgyro_catch_error('ERROR: Too many ions for NEO') 
  endif

  mu1 = sqrt(mi(1)/(me*loc_me_multiplier))

  ! Initialize NEO
  call neo_init(paths(i_r-1),gyro_comm)

  ! Simulation mode (dke solve vs. analytic)
  if (loc_neo_method == 1) then
     neo_sim_model_in = 0
  else
     neo_sim_model_in = 1
  end if

  ! Resolution
  if (loc_n_ion < 4) then
!     neo_n_energy_in = 10
     neo_n_energy_in = 5
  else
!     neo_n_energy_in = 6
     neo_n_energy_in = 5
  endif
  neo_n_xi_in = 11
  neo_n_theta_in = 11

  ! Geometry
  neo_equilibrium_model_in = 2
  neo_rmin_over_a_in       = r(i_r)/r_min
  neo_rmaj_over_a_in       = r_maj(i_r)/r_min
  neo_q_in                 = q(i_r)
  neo_shear_in             = s(i_r)
  neo_shift_in             = shift(i_r)
  neo_kappa_in             = kappa(i_r)
  neo_s_kappa_in           = s_kappa(i_r)
  neo_delta_in             = delta(i_r)
  neo_s_delta_in           = s_delta(i_r)

  neo_ipccw_in = tgyro_ipccw_in
  neo_btccw_in = tgyro_btccw_in

  neo_n_species_in = loc_n_ion+1
  neo_rho_star_in  = 0.001

  ! Electrons
  neo_z_2_in      = -1
  neo_mass_2_in   = 1.0/mu1**2
  neo_dens_2_in   = 1.0
  neo_temp_2_in   = 1.0/(ti(1,i_r)/te(i_r))
  neo_dlnndr_2_in = r_min*dlnnedr(i_r)
  neo_dlntdr_2_in = r_min*dlntedr(i_r)

  ! Main ions
  neo_z_1_in      = int(zi_vec(1))
  neo_mass_1_in   = 1.0
  neo_dens_1_in   = ni(1,i_r)/ne(i_r)
  neo_temp_1_in   = 1.0
  neo_dlnndr_1_in = r_min*dlnnidr(1,i_r)
  neo_dlntdr_1_in = r_min*dlntidr(1,i_r)
  neo_nu_1_in     = ((nue(i_r)*r_min/c_s(i_r))/&
       sqrt(ti(1,i_r)/te(i_r)))/(mu1*(ti(1,i_r)/te(i_r))**1.5)

  ! Ion #2
  if (neo_n_species_in > 2) then

     neo_z_3_in      = int(zi_vec(2))
     neo_mass_3_in   = 1.0/sqrt(mi(1)/mi(2))**2
     neo_dens_3_in   = ni(2,i_r)/ne(i_r)
     neo_dlnndr_3_in = r_min*dlnnidr(2,i_r)
     neo_dlntdr_3_in = r_min*dlntidr(2,i_r)

  endif

  ! Ion #3
  if (neo_n_species_in > 3) then

     neo_z_4_in      = int(zi_vec(3))
     neo_mass_4_in   = 1.0/sqrt(mi(1)/mi(3))**2
     neo_dens_4_in   = ni(3,i_r)/ne(i_r)
     neo_dlnndr_4_in = r_min*dlnnidr(3,i_r)
     neo_dlntdr_4_in = r_min*dlntidr(3,i_r)

  endif

  ! Ion #4
  if (neo_n_species_in > 4) then

     neo_z_5_in      = int(zi_vec(4))
     neo_mass_5_in   = 1.0/sqrt(mi(1)/mi(4))**2
     neo_dens_5_in   = ni(4,i_r)/ne(i_r)
     neo_dlnndr_5_in = r_min*dlnnidr(4,i_r)
     neo_dlntdr_5_in = r_min*dlntidr(4,i_r)

  endif

  ! Ion #5
  if (neo_n_species_in > 5) then

     neo_z_6_in      = int(zi_vec(5))
     neo_mass_6_in   = 1.0/sqrt(mi(1)/mi(5))**2
     neo_dens_6_in   = ni(5,i_r)/ne(i_r)
     neo_dlnndr_6_in = r_min*dlnnidr(5,i_r)
     neo_dlntdr_6_in = r_min*dlntidr(5,i_r)

  endif

  ! Rotation is always active
  neo_rotation_model_in = 2
  neo_omega_rot_in = u00(i_r) * r_min / r_maj(i_r) &
       / (c_s(i_r) * sqrt(ti(1,i_r)/te(i_r))) * (-1.0*tgyro_ipccw_in)
  neo_omega_rot_deriv_in = -gamma_p(i_r) * r_min**2 / r_maj(i_r) &
       / (c_s(i_r) * sqrt(ti(1,i_r)/te(i_r))) * (-1.0*tgyro_ipccw_in)

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
