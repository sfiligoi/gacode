subroutine neo_read_input

  use neo_globals

  implicit none

  integer :: is

  open(unit=1,file=trim(path)//'input.neo.gen',status='old')
  read(1,*) n_energy
  read(1,*) n_xi
  read(1,*) n_theta
  read(1,*) n_radial
  read(1,*) matsz_scalefac
  read(1,*) rmin_1_in
  read(1,*) rmin_2_in
  read(1,*) rmaj_in
  read(1,*) silent_flag
  read(1,*) sim_model
  read(1,*) equilibrium_model
  read(1,*) collision_model
  read(1,*) profile_model
  read(1,*) profile_erad0_model
  read(1,*) profile_equilibrium_model
  read(1,*) ipccw_in
  read(1,*) btccw_in
  read(1,*) te_ade_in
  read(1,*) ne_ade_in

  read(1,*) rotation_model
  read(1,*) omega_rot_in
  read(1,*) omega_rot_deriv_in

  read(1,*) spitzer_model
  read(1,*) epar0_spitzer

  read(1,*) n_species
 
  read(1,*) nu_1_in

  do is=1,6
     read(1,*) z_in(is)
     read(1,*) mass_in(is)
     read(1,*) dens_in(is)
     read(1,*) temp_in(is)
     read(1,*) dlnndr_in(is)
     read(1,*) dlntdr_in(is)
  enddo

  read(1,*) dphi0dr_in
  read(1,*) epar0_in
  read(1,*) q_in
  read(1,*) rho_in
  read(1,*) shat_in
  read(1,*) shift_in     
  read(1,*) kappa_in    
  read(1,*) s_kappa_in   
  read(1,*) delta_in        
  read(1,*) s_delta_in 
  read(1,*) zeta_in        
  read(1,*) s_zeta_in
  read(1,*) zmag_in        
  read(1,*) s_zmag_in
  read(1,*) profile_delta_scale
  read(1,*) profile_zeta_scale
  read(1,*) profile_zmag_scale

  read(1,*) profile_dlnndr_1_scale
  read(1,*) profile_dlnndr_2_scale
  read(1,*) profile_dlnndr_3_scale
  read(1,*) profile_dlnndr_4_scale
  read(1,*) profile_dlnndr_5_scale
  read(1,*) profile_dlnndr_6_scale

  read(1,*) profile_dlntdr_1_scale
  read(1,*) profile_dlntdr_2_scale
  read(1,*) profile_dlntdr_3_scale
  read(1,*) profile_dlntdr_4_scale
  read(1,*) profile_dlntdr_5_scale
  read(1,*) profile_dlntdr_6_scale


  close(1)

  ! GEO fourier coefficients are not yet available to read-in
  geo_ny_in = 0
  geo_yin_in(:,:) = 0.0

end subroutine neo_read_input
