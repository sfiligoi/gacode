!--------------------------------------------------------------
! tglf_dump_input.f90
!
! PURPOSE:
!  Complete dump of input parameters.
!
! NOTES:
!  Necessary for debugging exactly all input variables used in a simulation.
!
! 
!---------------------------------------------------------------

 subroutine tglf_dump_input

 use tglf_interface 
 implicit none
!
 integer :: funit=21
 integer :: gunit=22
 character(100) :: fname
 character(100) :: gname
!-------------------------
   fname=trim(tglf_path_in)//"tglf_dump_inputs.txt"
   open(unit=funit,file=trim(fname),status='replace')
   gname=trim(tglf_path_in)//"gen_input.dat"
   open(unit=gunit,file=trim(gname),status='replace')


 write(21,*) "tglf_use_transport_model_in = ",tglf_use_transport_model_in
 write(22,*) tglf_use_transport_model_in
 write(21,*) "tglf_geometry_flag_in = ",tglf_geometry_flag_in
 write(22,*) tglf_geometry_flag_in
 write(21,*) "tglf_write_wavefunction_flag_in = ",tglf_write_wavefunction_flag_in
 write(22,*) tglf_write_wavefunction_flag_in
 write(21,*) "tglf_sign_bt_in = ",tglf_sign_bt_in
 write(22,*) tglf_sign_bt_in
 write(21,*) "tglf_sign_it_in = ",tglf_sign_it_in
 write(22,*) tglf_sign_it_in
 write(21,*) "tglf_theta_trapped_in = ",tglf_theta_trapped_in
 write(22,*) tglf_theta_trapped_in
 write(21,*) "tglf_park_in = ",tglf_park_in
 write(22,*) tglf_park_in
 write(21,*) "tglf_ghat_in = ",tglf_ghat_in
 write(22,*) tglf_ghat_in
 write(21,*) "tglf_gchat_in = ",tglf_gchat_in
 write(22,*) tglf_gchat_in
 write(21,*) "tglf_wd_zero_in = ",tglf_wd_zero_in
 write(22,*) tglf_wd_zero_in
 write(21,*) "tglf_linsker_factor_in = ",tglf_linsker_factor_in
 write(22,*) tglf_linsker_factor_in
 write(21,*) "tglf_gradb_factor_in = ",tglf_gradb_factor_in
 write(22,*) tglf_gradb_factor_in
 write(21,*) "tglf_filter_in = ",tglf_filter_in
 write(22,*) tglf_filter_in
 write(21,*) "tglf_damp_psi_in = ",tglf_damp_psi_in
 write(22,*) tglf_damp_psi_in
 write(21,*) "tglf_damp_sig_in = ",tglf_damp_sig_in
 write(22,*) tglf_damp_sig_in
 write(21,*) "tglf_iflux_in = ",tglf_iflux_in
 write(22,*) tglf_iflux_in
 write(21,*) "tglf_use_bper_in = ",tglf_use_bper_in
 write(22,*) tglf_use_bper_in
 write(21,*) "tglf_use_bpar_in = ",tglf_use_bpar_in
 write(22,*) tglf_use_bpar_in
 write(21,*) "tglf_use_mhd_rule_in = ",tglf_use_mhd_rule_in
 write(22,*) tglf_use_mhd_rule_in
 write(21,*) "tglf_use_bisection_in = ",tglf_use_bisection_in
 write(22,*) tglf_use_bisection_in
 write(21,*) "tglf_ibranch_in = ",tglf_ibranch_in
 write(22,*) tglf_ibranch_in
 write(21,*) "tglf_nmodes_in = ",tglf_nmodes_in
 write(22,*) tglf_nmodes_in
 write(21,*) "tglf_nbasis_max_in = ",tglf_nbasis_max_in
 write(22,*) tglf_nbasis_max_in
 write(21,*) "tglf_nbasis_min_in = ",tglf_nbasis_min_in
 write(22,*) tglf_nbasis_min_in
 write(21,*) "tglf_nxgrid_in = ",tglf_nxgrid_in
 write(22,*) tglf_nxgrid_in
 write(21,*) "tglf_nky_in = ",tglf_nky_in
 write(22,*) tglf_nky_in
 write(21,*) "tglf_adiabatic_elec_in = ",tglf_adiabatic_elec_in
 write(22,*) tglf_adiabatic_elec_in
 write(21,*) "tglf_alpha_e_in = ",tglf_alpha_e_in
 write(22,*) tglf_alpha_e_in
 write(21,*) "tglf_alpha_p_in = ",tglf_alpha_p_in
 write(22,*) tglf_alpha_p_in
 write(21,*) "tglf_alpha_n_in = ",tglf_alpha_n_in
 write(22,*) tglf_alpha_n_in
 write(21,*) "tglf_alpha_t_in = ",tglf_alpha_t_in
 write(22,*) tglf_alpha_t_in
 write(21,*) "tglf_alpha_kx_e_in = ",tglf_alpha_kx_e_in
 write(22,*) tglf_alpha_kx_e_in
 write(21,*) "tglf_alpha_kx_p_in = ",tglf_alpha_kx_p_in
 write(22,*) tglf_alpha_kx_p_in
 write(21,*) "tglf_alpha_kx_n_in = ",tglf_alpha_kx_n_in
 write(22,*) tglf_alpha_kx_n_in
 write(21,*) "tglf_alpha_kx_t_in = ",tglf_alpha_kx_t_in
 write(22,*) tglf_alpha_kx_t_in
 write(21,*) "tglf_alpha_quench_in = ",tglf_alpha_quench_in
 write(22,*) tglf_alpha_quench_in
 write(21,*) "tglf_xnu_factor_in = ",tglf_xnu_factor_in
 write(22,*) tglf_xnu_factor_in
 write(21,*) "tglf_debye_factor_in = ",tglf_debye_factor_in
 write(22,*) tglf_debye_factor_in
 write(21,*) "tglf_etg_factor_in = ",tglf_etg_factor_in
 write(22,*) tglf_etg_factor_in
 write(21,*) "tglf_sat_rule_in = ",tglf_sat_rule_in
 write(22,*) tglf_sat_rule_in
 write(21,*) "tglf_kygrid_model_in = ",tglf_kygrid_model_in
 write(22,*) tglf_kygrid_model_in
 write(21,*) "tglf_xnu_model_in = ",tglf_xnu_model_in
 write(22,*) tglf_xnu_model_in
 write(21,*) "tglf_vpar_model_in = ",tglf_vpar_model_in
 write(22,*) tglf_vpar_model_in
 write(21,*) "tglf_vpar_shear_model_in = ",tglf_vpar_shear_model_in
 write(22,*) tglf_vpar_shear_model_in
 write(21,*) "tglf_ns_in = ",tglf_ns_in
 write(22,*) tglf_ns_in
 write(21,*) "tglf_mass_in(1) = ",tglf_mass_in(1)
 write(22,*) tglf_mass_in(1)
 write(21,*) "tglf_mass_in(2) = ",tglf_mass_in(2)
 write(22,*) tglf_mass_in(2)
 write(21,*) "tglf_mass_in(3) = ",tglf_mass_in(3)
 write(22,*) tglf_mass_in(3)
 write(21,*) "tglf_mass_in(4) = ",tglf_mass_in(4)
 write(22,*) tglf_mass_in(4)
 write(21,*) "tglf_mass_in(5) = ",tglf_mass_in(5)
 write(22,*) tglf_mass_in(5)
 write(21,*) "tglf_mass_in(6) = ",tglf_mass_in(6)
 write(22,*) tglf_mass_in(6)
 write(21,*) "tglf_zs_in(1) = ",tglf_zs_in(1)
 write(22,*) tglf_zs_in(1)
 write(21,*) "tglf_zs_in(2) = ",tglf_zs_in(2)
 write(22,*) tglf_zs_in(2)
 write(21,*) "tglf_zs_in(3) = ",tglf_zs_in(3)
 write(22,*) tglf_zs_in(3)
 write(21,*) "tglf_zs_in(4) = ",tglf_zs_in(4)
 write(22,*) tglf_zs_in(4)
 write(21,*) "tglf_zs_in(5) = ",tglf_zs_in(5)
 write(22,*) tglf_zs_in(5)
 write(21,*) "tglf_zs_in(6) = ",tglf_zs_in(6)
 write(22,*) tglf_zs_in(6)
 write(21,*) "tglf_ky_in = ",tglf_ky_in
 write(22,*) tglf_ky_in
 write(21,*) "tglf_width_in = ",tglf_width_in
 write(22,*) tglf_width_in
 write(21,*) "tglf_width_min_in = ",tglf_width_min_in
 write(22,*) tglf_width_min_in
 write(21,*) "tglf_nwidth_in = ",tglf_nwidth_in
 write(22,*) tglf_nwidth_in
 write(21,*) "tglf_find_width_in = ",tglf_find_width_in
 write(22,*) tglf_find_width_in
 write(21,*) "tglf_rlns_in(1) = ",tglf_rlns_in(1)
 write(22,*) tglf_rlns_in(1)
 write(21,*) "tglf_rlns_in(2) = ",tglf_rlns_in(2)
 write(22,*) tglf_rlns_in(2)
 write(21,*) "tglf_rlns_in(3) = ",tglf_rlns_in(3)
 write(22,*) tglf_rlns_in(3)
 write(21,*) "tglf_rlns_in(4) = ",tglf_rlns_in(4)
 write(22,*) tglf_rlns_in(4)
 write(21,*) "tglf_rlns_in(5) = ",tglf_rlns_in(5)
 write(22,*) tglf_rlns_in(5)
 write(21,*) "tglf_rlns_in(6) = ",tglf_rlns_in(6)
 write(22,*) tglf_rlns_in(6)
 write(21,*) "tglf_rlts_in(1) = ",tglf_rlts_in(1)
 write(22,*) tglf_rlts_in(1)
 write(21,*) "tglf_rlts_in(2) = ",tglf_rlts_in(2)
 write(22,*) tglf_rlts_in(2)
 write(21,*) "tglf_rlts_in(3) = ",tglf_rlts_in(3)
 write(22,*) tglf_rlts_in(3)
 write(21,*) "tglf_rlts_in(4) = ",tglf_rlts_in(4)
 write(22,*) tglf_rlts_in(4)
 write(21,*) "tglf_rlts_in(5) = ",tglf_rlts_in(5)
 write(22,*) tglf_rlts_in(5)
 write(21,*) "tglf_rlts_in(6) = ",tglf_rlts_in(6)
 write(22,*) tglf_rlts_in(6)
 write(21,*) "tglf_vpar_shear_in(1) = ",tglf_vpar_shear_in(1)
 write(22,*) tglf_vpar_shear_in(1)
 write(21,*) "tglf_vpar_shear_in(2) = ",tglf_vpar_shear_in(2)
 write(22,*) tglf_vpar_shear_in(2)
 write(21,*) "tglf_vpar_shear_in(3) = ",tglf_vpar_shear_in(3)
 write(22,*) tglf_vpar_shear_in(3)
 write(21,*) "tglf_vpar_shear_in(4) = ",tglf_vpar_shear_in(4)
 write(22,*) tglf_vpar_shear_in(4)
 write(21,*) "tglf_vpar_shear_in(5) = ",tglf_vpar_shear_in(5)
 write(22,*) tglf_vpar_shear_in(5)
 write(21,*) "tglf_vpar_shear_in(6) = ",tglf_vpar_shear_in(6)
 write(22,*) tglf_vpar_shear_in(6)
 write(21,*) "tglf_vexb_shear_in = ",tglf_vexb_shear_in
 write(22,*) tglf_vexb_shear_in
 write(21,*) "tglf_vns_shear_in(1) = ",tglf_vns_shear_in(1)
 write(22,*) tglf_vns_shear_in(1)
 write(21,*) "tglf_vns_shear_in(2) = ",tglf_vns_shear_in(2)
 write(22,*) tglf_vns_shear_in(2)
 write(21,*) "tglf_vns_shear_in(3) = ",tglf_vns_shear_in(3)
 write(22,*) tglf_vns_shear_in(3)
 write(21,*) "tglf_vns_shear_in(4) = ",tglf_vns_shear_in(4)
 write(22,*) tglf_vns_shear_in(4)
 write(21,*) "tglf_vns_shear_in(5) = ",tglf_vns_shear_in(5)
 write(22,*) tglf_vns_shear_in(5)
 write(21,*) "tglf_vns_shear_in(6) = ",tglf_vns_shear_in(6)
 write(22,*) tglf_vns_shear_in(6)
 write(21,*) "tglf_vts_shear_in(1) = ",tglf_vts_shear_in(1)
 write(22,*) tglf_vts_shear_in(1)
 write(21,*) "tglf_vts_shear_in(2) = ",tglf_vts_shear_in(2)
 write(22,*) tglf_vts_shear_in(2)
 write(21,*) "tglf_vts_shear_in(3) = ",tglf_vts_shear_in(3)
 write(22,*) tglf_vts_shear_in(3)
 write(21,*) "tglf_vts_shear_in(4) = ",tglf_vts_shear_in(4)
 write(22,*) tglf_vts_shear_in(4)
 write(21,*) "tglf_vts_shear_in(5) = ",tglf_vts_shear_in(5)
 write(22,*) tglf_vts_shear_in(5)
 write(21,*) "tglf_vts_shear_in(6) = ",tglf_vts_shear_in(6)
 write(22,*) tglf_vts_shear_in(6)
 write(21,*) "tglf_taus_in(1) = ",tglf_taus_in(1)
 write(22,*) tglf_taus_in(1)
 write(21,*) "tglf_taus_in(2) = ",tglf_taus_in(2)
 write(22,*) tglf_taus_in(2)
 write(21,*) "tglf_taus_in(3) = ",tglf_taus_in(3)
 write(22,*) tglf_taus_in(3)
 write(21,*) "tglf_taus_in(4) = ",tglf_taus_in(4)
 write(22,*) tglf_taus_in(4)
 write(21,*) "tglf_taus_in(5) = ",tglf_taus_in(5)
 write(22,*) tglf_taus_in(5)
 write(21,*) "tglf_taus_in(6) = ",tglf_taus_in(6)
 write(22,*) tglf_taus_in(6)
 write(21,*) "tglf_as_in(1) = ",tglf_as_in(1)
 write(22,*) tglf_as_in(1)
 write(21,*) "tglf_as_in(2) = ",tglf_as_in(2)
 write(22,*) tglf_as_in(2)
 write(21,*) "tglf_as_in(3) = ",tglf_as_in(3)
 write(22,*) tglf_as_in(3)
 write(21,*) "tglf_as_in(4) = ",tglf_as_in(4)
 write(22,*) tglf_as_in(4)
 write(21,*) "tglf_as_in(5) = ",tglf_as_in(5)
 write(22,*) tglf_as_in(5)
 write(21,*) "tglf_as_in(6) = ",tglf_as_in(6)
 write(22,*) tglf_as_in(6)
 write(21,*) "tglf_vpar_in(1) = ",tglf_vpar_in(1)
 write(22,*) tglf_vpar_in(1)
 write(21,*) "tglf_vpar_in(2) = ",tglf_vpar_in(2)
 write(22,*) tglf_vpar_in(2)
 write(21,*) "tglf_vpar_in(3) = ",tglf_vpar_in(3)
 write(22,*) tglf_vpar_in(3)
 write(21,*) "tglf_vpar_in(4) = ",tglf_vpar_in(4)
 write(22,*) tglf_vpar_in(4)
 write(21,*) "tglf_vpar_in(5) = ",tglf_vpar_in(5)
 write(22,*) tglf_vpar_in(5)
 write(21,*) "tglf_vpar_in(6) = ",tglf_vpar_in(6)
 write(22,*) tglf_vpar_in(6)
 write(21,*) "tglf_betae_in = ",tglf_betae_in
 write(22,*) tglf_betae_in
 write(21,*) "tglf_xnue_in = ",tglf_xnue_in
 write(22,*) tglf_xnue_in
 write(21,*) "tglf_zeff_in = ",tglf_zeff_in
 write(22,*) tglf_zeff_in
 write(21,*) "tglf_debye_in = ",tglf_debye_in
 write(22,*) tglf_debye_in
 write(21,*) "tglf_new_eikonal_in = ",tglf_new_eikonal_in
 write(22,*) tglf_new_eikonal_in
 write(21,*) "tglf_rmin_sa_in = ",tglf_rmin_sa_in
 write(22,*) tglf_rmin_sa_in
 write(21,*) "tglf_rmaj_sa_in = ",tglf_rmaj_sa_in
 write(22,*) tglf_rmaj_sa_in
 write(21,*) "tglf_q_sa_in = ",tglf_q_sa_in
 write(22,*) tglf_q_sa_in
 write(21,*) "tglf_shat_sa_in = ",tglf_shat_sa_in
 write(22,*) tglf_shat_sa_in
 write(21,*) "tglf_alpha_sa_in = ",tglf_alpha_sa_in
 write(22,*) tglf_alpha_sa_in
 write(21,*) "tglf_xwell_sa_in = ",tglf_xwell_sa_in
 write(22,*) tglf_xwell_sa_in
 write(21,*) "tglf_theta0_sa_in = ",tglf_theta0_sa_in
 write(22,*) tglf_theta0_sa_in
 write(21,*) "tglf_b_model_sa_in = ",tglf_b_model_sa_in
 write(22,*) tglf_b_model_sa_in
 write(21,*) "tglf_ft_model_sa_in = ",tglf_ft_model_sa_in
 write(22,*) tglf_ft_model_sa_in
 write(21,*) "tglf_rmin_loc_in = ",tglf_rmin_loc_in
 write(22,*) tglf_rmin_loc_in
 write(21,*) "tglf_rmaj_loc_in = ",tglf_rmaj_loc_in
 write(22,*) tglf_rmaj_loc_in
 write(21,*) "tglf_zmaj_loc_in = ",tglf_zmaj_loc_in
 write(22,*) tglf_zmaj_loc_in
 write(21,*) "tglf_drmindx_loc_in = ",tglf_drmindx_loc_in
 write(22,*) tglf_drmindx_loc_in
 write(21,*) "tglf_drmajdx_loc_in = ",tglf_drmajdx_loc_in
 write(22,*) tglf_drmajdx_loc_in
 write(21,*) "tglf_dzmajdx_loc_in = ",tglf_dzmajdx_loc_in
 write(22,*) tglf_dzmajdx_loc_in
 write(21,*) "tglf_kappa_loc_in = ",tglf_kappa_loc_in
 write(22,*) tglf_kappa_loc_in
 write(21,*) "tglf_s_kappa_loc_in = ",tglf_s_kappa_loc_in
 write(22,*) tglf_s_kappa_loc_in
 write(21,*) "tglf_delta_loc_in = ",tglf_delta_loc_in
 write(22,*) tglf_delta_loc_in
 write(21,*) "tglf_s_delta_loc_in = ",tglf_s_delta_loc_in
 write(22,*) tglf_s_delta_loc_in
 write(21,*) "tglf_zeta_loc_in = ",tglf_zeta_loc_in
 write(22,*) tglf_zeta_loc_in
 write(21,*) "tglf_s_zeta_loc_in = ",tglf_s_zeta_loc_in
 write(22,*) tglf_s_zeta_loc_in
 write(21,*) "tglf_q_loc_in = ",tglf_q_loc_in
 write(22,*) tglf_q_loc_in
 write(21,*) "tglf_q_prime_loc_in = ",tglf_q_prime_loc_in
 write(22,*) tglf_q_prime_loc_in
 write(21,*) "tglf_p_prime_loc_in = ",tglf_p_prime_loc_in
 write(22,*) tglf_p_prime_loc_in
  close(funit)
!
 end subroutine tglf_dump_input
!--------------------------------------------------------

