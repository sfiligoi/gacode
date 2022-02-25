!--------------------------------------------------------------
! gftm_read_input.f90
!
! PURPOSE:
!  Complete read of gftm input parameters as specified by
!  gftm_interface (but ignoring Fourier and ELITE geometry
!  coefficients).
!--------------------------------------------------------------

subroutine gftm_read_input

  use gftm_interface

  implicit none

  open(unit=1,file=trim(gftm_path_in)//'input.gftm.gen',status='old')

  read(1,*) gftm_units_in
  read(1,*) gftm_use_transport_model_in
  read(1,*) gftm_geometry_flag_in
  read(1,*) gftm_write_wavefunction_flag_in

  ! Data passed to: put_signs
  read(1,*) gftm_sign_bt_in
  read(1,*) gftm_sign_it_in

  ! Data passed to: put_rare_switches
  read(1,*) gftm_theta_trapped_in
  read(1,*) gftm_wdia_trapped_in
  read(1,*) gftm_park_in
  read(1,*) gftm_ghat_in
  read(1,*) gftm_gchat_in
  read(1,*) gftm_wd_zero_in
  read(1,*) gftm_linsker_factor_in
  read(1,*) gftm_gradb_factor_in
  read(1,*) gftm_filter_in
  read(1,*) gftm_damp_psi_in
  read(1,*) gftm_damp_sig_in

  ! Data passed to: put_switches
  read(1,*) gftm_iflux_in
  read(1,*) gftm_use_bper_in
  read(1,*) gftm_use_bpar_in
  read(1,*) gftm_use_mhd_rule_in
  read(1,*) gftm_use_bisection_in
  read(1,*) gftm_use_inboard_detrapped_in
  read(1,*) gftm_ibranch_in
  read(1,*) gftm_nmodes_in
  read(1,*) gftm_nbasis_max_in
  read(1,*) gftm_nbasis_min_in
  read(1,*) gftm_nxgrid_in
  read(1,*) gftm_nky_in
  read(1,*) gftm_use_ave_ion_grid_in

  ! Data passed to: put_model_parameters
  read(1,*) gftm_adiabatic_elec_in
  read(1,*) gftm_alpha_mach_in
  read(1,*) gftm_alpha_e_in
  read(1,*) gftm_alpha_p_in
  read(1,*) gftm_alpha_quench_in
  read(1,*) gftm_alpha_zf_in
  read(1,*) gftm_xnu_factor_in
  read(1,*) gftm_debye_factor_in
  read(1,*) gftm_etg_factor_in
  read(1,*) gftm_rlnp_cutoff_in
  read(1,*) gftm_sat_rule_in
  read(1,*) gftm_kygrid_model_in
  read(1,*) gftm_xnu_model_in
  read(1,*) gftm_vpar_model_in
  read(1,*) gftm_vpar_shear_model_in

  ! Data passed to: put_species
  read(1,*) gftm_ns_in
  read(1,*) gftm_mass_in(1)
  read(1,*) gftm_mass_in(2)
  read(1,*) gftm_mass_in(3)
  read(1,*) gftm_mass_in(4)
  read(1,*) gftm_mass_in(5)
  read(1,*) gftm_mass_in(6)
  read(1,*) gftm_mass_in(7)
  read(1,*) gftm_zs_in(1)
  read(1,*) gftm_zs_in(2)
  read(1,*) gftm_zs_in(3)
  read(1,*) gftm_zs_in(4)
  read(1,*) gftm_zs_in(5)
  read(1,*) gftm_zs_in(6)
  read(1,*) gftm_zs_in(7)

  ! Data passed to: put_kys
  read(1,*) gftm_ky_in

  ! Data passed to: put_gaussian_width
  read(1,*) gftm_width_in
  read(1,*) gftm_width_min_in
  read(1,*) gftm_nwidth_in
  read(1,*) gftm_find_width_in

  ! Data passed to: put_gradients
  read(1,*) gftm_rlns_in(1)
  read(1,*) gftm_rlns_in(2)
  read(1,*) gftm_rlns_in(3)
  read(1,*) gftm_rlns_in(4)
  read(1,*) gftm_rlns_in(5)
  read(1,*) gftm_rlns_in(6)
  read(1,*) gftm_rlns_in(7)
  read(1,*) gftm_rlts_in(1)
  read(1,*) gftm_rlts_in(2)
  read(1,*) gftm_rlts_in(3)
  read(1,*) gftm_rlts_in(4)
  read(1,*) gftm_rlts_in(5)
  read(1,*) gftm_rlts_in(6)
  read(1,*) gftm_rlts_in(7)
  read(1,*) gftm_vpar_shear_in(1)
  read(1,*) gftm_vpar_shear_in(2)
  read(1,*) gftm_vpar_shear_in(3)
  read(1,*) gftm_vpar_shear_in(4)
  read(1,*) gftm_vpar_shear_in(5)
  read(1,*) gftm_vpar_shear_in(6)
  read(1,*) gftm_vpar_shear_in(7)
  read(1,*) gftm_vexb_shear_in

  ! Data passed to: put_profile_shear
  read(1,*) gftm_vns_shear_in(1)
  read(1,*) gftm_vns_shear_in(2)
  read(1,*) gftm_vns_shear_in(3)
  read(1,*) gftm_vns_shear_in(4)
  read(1,*) gftm_vns_shear_in(5)
  read(1,*) gftm_vns_shear_in(6)
  read(1,*) gftm_vns_shear_in(7)
  read(1,*) gftm_vts_shear_in(1)
  read(1,*) gftm_vts_shear_in(2)
  read(1,*) gftm_vts_shear_in(3)
  read(1,*) gftm_vts_shear_in(4)
  read(1,*) gftm_vts_shear_in(5)
  read(1,*) gftm_vts_shear_in(6)
  read(1,*) gftm_vts_shear_in(7)

  ! Data passed to: put_averages
  read(1,*) gftm_taus_in(1)
  read(1,*) gftm_taus_in(2)
  read(1,*) gftm_taus_in(3)
  read(1,*) gftm_taus_in(4)
  read(1,*) gftm_taus_in(5)
  read(1,*) gftm_taus_in(6)
  read(1,*) gftm_taus_in(7)
  read(1,*) gftm_as_in(1)
  read(1,*) gftm_as_in(2)
  read(1,*) gftm_as_in(3)
  read(1,*) gftm_as_in(4)
  read(1,*) gftm_as_in(5)
  read(1,*) gftm_as_in(6)
  read(1,*) gftm_as_in(7)
  read(1,*) gftm_vpar_in(1)
  read(1,*) gftm_vpar_in(2)
  read(1,*) gftm_vpar_in(3)
  read(1,*) gftm_vpar_in(4)
  read(1,*) gftm_vpar_in(5)
  read(1,*) gftm_vpar_in(6)
  read(1,*) gftm_vpar_in(7)
  read(1,*) gftm_vexb_in
  read(1,*) gftm_betae_in
  read(1,*) gftm_xnue_in
  read(1,*) gftm_zeff_in
  read(1,*) gftm_debye_in

  ! Data passed to: put_eikonal
  read(1,*) gftm_new_eikonal_in

  ! Data passed to: put_s_alpha_geometry
  read(1,*) gftm_rmin_sa_in
  read(1,*) gftm_rmaj_sa_in
  read(1,*) gftm_q_sa_in
  read(1,*) gftm_shat_sa_in
  read(1,*) gftm_alpha_sa_in
  read(1,*) gftm_xwell_sa_in
  read(1,*) gftm_theta0_sa_in
  read(1,*) gftm_b_model_sa_in
  read(1,*) gftm_ft_model_sa_in

  ! Data passed to: put_Miller_geometry
  read(1,*) gftm_rmin_loc_in
  read(1,*) gftm_rmaj_loc_in
  read(1,*) gftm_zmaj_loc_in
  read(1,*) gftm_drmindx_loc_in
  read(1,*) gftm_drmajdx_loc_in
  read(1,*) gftm_dzmajdx_loc_in
  read(1,*) gftm_kappa_loc_in
  read(1,*) gftm_s_kappa_loc_in
  read(1,*) gftm_delta_loc_in
  read(1,*) gftm_s_delta_loc_in
  read(1,*) gftm_zeta_loc_in
  read(1,*) gftm_s_zeta_loc_in
  read(1,*) gftm_q_loc_in
  read(1,*) gftm_q_prime_loc_in
  read(1,*) gftm_p_prime_loc_in
  read(1,*) gftm_kx0_loc_in

  ! Set threshold for gftm-NN execution versus full gftm calculation
  read(1,*) gftm_nn_max_error_in

  close(1)

end subroutine gftm_read_input
