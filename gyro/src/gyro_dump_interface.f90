!--------------------------------------------------------------
! gyro_dump_interface.f90
!
! PURPOSE:
!  Complete dump of interface parameters matched with.
!
! NOTES:
!  Also necessary for debugging exactly all interface variables used in a simulation.
!
! 
!---------------------------------------------------------------

subroutine gyro_dump_interface

  use gyro_globals
  use gyro_interface
  implicit none
  !
  integer :: funit=21
  character(100) :: fname
  !-------------------------
  if (i_proc == 0) THEN
     fname=trim(path)//'out.gyro.interfacedump'
     open(unit=funit,file=trim(fname),status='replace')


     call dumpIntInterface(21,"    gyro_radial_grid_in ",    gyro_radial_grid_in , &
          "n_x",n_x)
     call dumpIntInterface(21,"    gyro_orbit_grid_in ",    gyro_orbit_grid_in , &
          "n_theta_section",n_theta_section)
     call dumpIntInterface(21,"    gyro_pass_grid_in ",    gyro_pass_grid_in , &
          "n_pass",n_pass)
     call dumpIntInterface(21,"    gyro_trap_grid_in ",    gyro_trap_grid_in , &
          "n_trap",n_trap)
     call dumpIntInterface(21,"    gyro_energy_grid_in ",    gyro_energy_grid_in , &
          "n_energy",n_energy)
     call dumpIntInterface(21,"    gyro_variable_egrid_flag_in ",    gyro_variable_egrid_flag_in , &
          "variable_egrid_flag",variable_egrid_flag)
     call dumpIntInterface(21,"    gyro_blend_grid_in ",    gyro_blend_grid_in , &
          "n_blend",n_blend)
     call dumpIntInterface(21,"    gyro_blend_fit_order_in ",    gyro_blend_fit_order_in , &
          "blend_fit_order",blend_fit_order)
     call dumpIntInterface(21,"    gyro_theta_plot_in ",    gyro_theta_plot_in , &
          "n_theta_plot",n_theta_plot)
     call dumpIntInterface(21,"    gyro_toroidal_min_in ",    gyro_toroidal_min_in , &
          "n0",n0)
     call dumpIntInterface(21,"    gyro_toroidal_grid_in ",    gyro_toroidal_grid_in , &
          "n_n",n_n)
     call dumpIntInterface(21,"    gyro_toroidal_sep_in ",    gyro_toroidal_sep_in , &
          "d_n",d_n)
     call dumpIntInterface(21,"    gyro_toroidal_reference_in ",    gyro_toroidal_reference_in , &
          "n_ref",n_ref)
     call dumpRealInterface(21,"    gyro_safety_factor_in ",    gyro_safety_factor_in , &
          "q0",q0)
     call dumpRealInterface(21,"    gyro_shear_in ",    gyro_shear_in , &
          "s0",s0)
     call dumpRealInterface(21,"    gyro_delta_in ",    gyro_delta_in , &
          "delta0",delta0)
     call dumpRealInterface(21,"    gyro_zeta_in ",    gyro_zeta_in , &
          "zeta0",zeta0)
     call dumpRealInterface(21,"    gyro_s_delta_in ",    gyro_s_delta_in , &
          "s_delta0",s_delta0)
     call dumpRealInterface(21,"    gyro_s_zeta_in ",    gyro_s_zeta_in , &
          "s_zeta0",s_zeta0)
     call dumpRealInterface(21,"    gyro_kappa_in ",    gyro_kappa_in , &
          "kappa0",kappa0)
     call dumpRealInterface(21,"    gyro_s_kappa_in ",    gyro_s_kappa_in , &
          "s_kappa0",s_kappa0)
     call dumpRealInterface(21,"    gyro_shift_in ",    gyro_shift_in , &
          "drmaj0",drmaj0)
     call dumpRealInterface(21,"    gyro_aspect_ratio_in ",    gyro_aspect_ratio_in , &
          "r_maj",r_maj)
     call dumpRealInterface(21,"    gyro_radius_in ",    gyro_radius_in , &
          "r0",r0)
     call dumpRealInterface(21,"    gyro_box_multiplier_in ",    gyro_box_multiplier_in , &
          "box_multiplier",box_multiplier)
     call dumpRealInterface(21,"    gyro_rho_star_in ",    gyro_rho_star_in , &
          "rho_star",rho_star)
     call dumpRealInterface(21,"    gyro_time_step_in ",    gyro_time_step_in , &
          "dt",dt)
     call dumpRealInterface(21,"    gyro_time_max_in ",    gyro_time_max_in , &
          "time_max",time_max)
     call dumpIntInterface(21,"    gyro_time_skip_in ",    gyro_time_skip_in , &
          "time_skip",time_skip)
     call dumpIntInterface(21,"    gyro_boundary_method_in ",    gyro_boundary_method_in , &
          "boundary_method",boundary_method)
     call dumpIntInterface(21,"    gyro_radial_derivative_band_in ",    gyro_radial_derivative_band_in , &
          "m_dx",m_dx)
     call dumpIntInterface(21,"    gyro_radial_gyro_band_in ",    gyro_radial_gyro_band_in , &
          "m_gyro",m_gyro)
     call dumpIntInterface(21,"    gyro_explicit_damp_grid_in ",    gyro_explicit_damp_grid_in , &
          "n_explicit_damp",n_explicit_damp)
     call dumpRealInterface(21,"    gyro_explicit_damp_in ",    gyro_explicit_damp_in , &
          "explicit_damp",explicit_damp)
     call dumpRealInterface(21,"    gyro_explicit_damp_elec_in ",    gyro_explicit_damp_elec_in , &
          "explicit_damp_elec",explicit_damp_elec)
     call dumpIntInterface(21,"    gyro_debug_flag_in ",    gyro_debug_flag_in , &
          "debug_flag",debug_flag)
     call dumpIntInterface(21,"    gyro_nonlinear_flag_in ",    gyro_nonlinear_flag_in , &
          "nonlinear_flag",nonlinear_flag)
     call dumpRealInterface(21,"    gyro_amp_phi_n_in ",    gyro_amp_phi_n_in , &
          "amp_n",amp_n)
     call dumpRealInterface(21,"    gyro_amp_phi_0_in ",    gyro_amp_phi_0_in , &
          "amp_0",amp_0)
     call dumpRealInterface(21,"    gyro_radial_upwind_in ",    gyro_radial_upwind_in , &
          "radial_upwind",radial_upwind)
     call dumpIntInterface(21,"    gyro_electron_method_in ",    gyro_electron_method_in , &
          "electron_method",electron_method)
     call dumpIntInterface(21,"    gyro_radial_profile_method_in ",    gyro_radial_profile_method_in , &
          "radial_profile_method",radial_profile_method)
     call dumpIntInterface(21,"    gyro_plot_u_flag_in ",    gyro_plot_u_flag_in , &
          "plot_u_flag",plot_u_flag)
     call dumpIntInterface(21,"    gyro_plot_epar_flag_in ",    gyro_plot_epar_flag_in , &
          "plot_epar_flag",plot_epar_flag)
     call dumpIntInterface(21,"    gyro_plot_n_flag_in ",    gyro_plot_n_flag_in , &
          "plot_n_flag",plot_n_flag)
     call dumpIntInterface(21,"    gyro_plot_e_flag_in ",    gyro_plot_e_flag_in , &
          "plot_e_flag",plot_e_flag)
     call dumpIntInterface(21,"    gyro_plot_v_flag_in ",    gyro_plot_v_flag_in , &
          "plot_v_flag",plot_v_flag)
     call dumpRealInterface(21,"    gyro_z_eff_in ",    gyro_z_eff_in , &
          "z_eff",z_eff)
     call dumpRealInterface(21,"    gyro_betae_unit_in ",    gyro_betae_unit_in , &
          "betae_unit",betae_unit)
     call dumpRealInterface(21,"    gyro_ampere_scale_in ",    gyro_ampere_scale_in , &
          "ampere_scale",ampere_scale)
     call dumpIntInterface(21,"    gyro_n_field_in ",    gyro_n_field_in , &
          "n_field",n_field)
     call dumpIntInterface(21,"    gyro_source_flag_in ",    gyro_source_flag_in , &
          "source_flag",source_flag)
     call dumpRealInterface(21,"    gyro_nu_source_in ",    gyro_nu_source_in , &
          "nu_source",nu_source)
     call dumpIntInterface(21,"    gyro_verbose_flag_in ",    gyro_verbose_flag_in , &
          "verbose_flag",verbose_flag)
     call dumpIntInterface(21,"    gyro_nonuniform_grid_flag_in ",    gyro_nonuniform_grid_flag_in , &
          "nonuniform_grid_flag",nonuniform_grid_flag)
     call dumpRealInterface(21,"    gyro_s_grid_in ",    gyro_s_grid_in , &
          "s_grid",s_grid)
     call dumpIntInterface(21,"    gyro_rotation_theory_method_in ",    gyro_rotation_theory_method_in , &
          "rotation_theory_method",rotation_theory_method)
     call dumpRealInterface(21,"    gyro_gamma_e_in ",    gyro_gamma_e_in , &
          "gamma_e",gamma_e)
     call dumpRealInterface(21,"    gyro_nu_ei_in ",    gyro_nu_ei_in , &
          "nu_ei",nu_ei)
     call dumpRealInterface(21,"    gyro_nu_ei_scale_in ",    gyro_nu_ei_scale_in , &
          "nu_ei_scale",nu_ei_scale)
     call dumpRealInterface(21,"    gyro_nu_ii_scale_in ",    gyro_nu_ii_scale_in , &
          "nu_ii_scale",nu_ii_scale)
     call dumpRealInterface(21,"    gyro_nu_i_krook_in ",    gyro_nu_i_krook_in , &
          "nu_i_krook",nu_i_krook)
     call dumpRealInterface(21,"    gyro_plot_filter_in ",    gyro_plot_filter_in , &
          "plot_filter",plot_filter)
     call dumpIntInterface(21,"    gyro_n_source_in ",    gyro_n_source_in , &
          "n_source",n_source)
     call dumpIntInterface(21,"    gyro_flat_profile_flag_in ",    gyro_flat_profile_flag_in , &
          "flat_profile_flag",flat_profile_flag)
     call dumpIntInterface(21,"    gyro_density_method_in ",    gyro_density_method_in , &
          "density_method",density_method)
     call dumpIntInterface(21,"    gyro_integrator_method_in ",    gyro_integrator_method_in , &
          "integrator_method",integrator_method)
     call dumpRealInterface(21,"    gyro_energy_max_1_in ",    gyro_energy_max_1_in , &
          "energy_max_vec(1)",energy_max_vec(1))
     call dumpRealInterface(21,"    gyro_energy_max_2_in ",    gyro_energy_max_2_in , &
          "energy_max_vec(2)",energy_max_vec(2))
     call dumpRealInterface(21,"    gyro_energy_max_3_in ",    gyro_energy_max_3_in , &
          "energy_max_vec(3)",energy_max_vec(3))
     call dumpRealInterface(21,"    gyro_energy_max_4_in ",    gyro_energy_max_4_in , &
          "energy_max_vec(4)",energy_max_vec(4))
     call dumpRealInterface(21,"    gyro_energy_max_5_in ",    gyro_energy_max_5_in , &
          "energy_max_vec(5)",energy_max_vec(5))
     call dumpRealInterface(21,"    gyro_energy_max_electron_in ",    gyro_energy_max_electron_in , &
          "energy_max_vec(0)",energy_max_vec(0))
     call dumpRealInterface(21,"    gyro_mu_1_in ",    gyro_mu_1_in , &
          "mu_vec(1)",mu_vec(1))
     call dumpRealInterface(21,"    gyro_mu_2_in ",    gyro_mu_2_in , &
          "mu_vec(2)",mu_vec(2))
     call dumpRealInterface(21,"    gyro_mu_3_in ",    gyro_mu_3_in , &
          "mu_vec(3)",mu_vec(3))
     call dumpRealInterface(21,"    gyro_mu_4_in ",    gyro_mu_4_in , &
          "mu_vec(4)",mu_vec(4))
     call dumpRealInterface(21,"    gyro_mu_5_in ",    gyro_mu_5_in , &
          "mu_vec(5)",mu_vec(5))
     call dumpRealInterface(21,"    gyro_mu_electron_in ",    gyro_mu_electron_in , &
          "mu_vec(0)",mu_vec(0))
     call dumpRealInterface(21,"    gyro_dlnndr_in ",    gyro_dlnndr_in , &
          "dlnndr_vec(1)",dlnndr_vec(1))
     call dumpRealInterface(21,"    gyro_dlnndr_2_in ",    gyro_dlnndr_2_in , &
          "dlnndr_vec(2)",dlnndr_vec(2))
     call dumpRealInterface(21,"    gyro_dlnndr_3_in ",    gyro_dlnndr_3_in , &
          "dlnndr_vec(3)",dlnndr_vec(3))
     call dumpRealInterface(21,"    gyro_dlnndr_4_in ",    gyro_dlnndr_4_in , &
          "dlnndr_vec(4)",dlnndr_vec(4))
     call dumpRealInterface(21,"    gyro_dlnndr_5_in ",    gyro_dlnndr_5_in , &
          "dlnndr_vec(5)",dlnndr_vec(5))
     call dumpRealInterface(21,"    gyro_dlnndr_electron_in ",    gyro_dlnndr_electron_in , &
          "dlnndr_vec(0)",dlnndr_vec(0))
     call dumpRealInterface(21,"    gyro_dlntdr_in ",    gyro_dlntdr_in , &
          "dlntdr_vec(1)",dlntdr_vec(1))
     call dumpRealInterface(21,"    gyro_dlntdr_2_in ",    gyro_dlntdr_2_in , &
          "dlntdr_vec(2)",dlntdr_vec(2))
     call dumpRealInterface(21,"    gyro_dlntdr_3_in ",    gyro_dlntdr_3_in , &
          "dlntdr_vec(3)",dlntdr_vec(3))
     call dumpRealInterface(21,"    gyro_dlntdr_4_in ",    gyro_dlntdr_4_in , &
          "dlntdr_vec(4)",dlntdr_vec(4))
     call dumpRealInterface(21,"    gyro_dlntdr_5_in ",    gyro_dlntdr_5_in , &
          "dlntdr_vec(5)",dlntdr_vec(5))
     call dumpRealInterface(21,"    gyro_dlntdr_electron_in ",    gyro_dlntdr_electron_in , &
          "dlntdr_vec(0)",dlntdr_vec(0))
     call dumpRealInterface(21,"    gyro_ni_over_ne_in ",    gyro_ni_over_ne_in , &
          "n_vec(1)",n_vec(1))
     call dumpRealInterface(21,"    gyro_ni_over_ne_2_in ",    gyro_ni_over_ne_2_in , &
          "n_vec(2)",n_vec(2))
     call dumpRealInterface(21,"    gyro_ni_over_ne_3_in ",    gyro_ni_over_ne_3_in , &
          "n_vec(3)",n_vec(3))
     call dumpRealInterface(21,"    gyro_ni_over_ne_4_in ",    gyro_ni_over_ne_4_in , &
          "n_vec(4)",n_vec(4))
     call dumpRealInterface(21,"    gyro_ni_over_ne_5_in ",    gyro_ni_over_ne_5_in , &
          "n_vec(5)",n_vec(5))
     call dumpRealInterface(21,"    gyro_ti_over_te_in ",    gyro_ti_over_te_in , &
          "t_vec(1)",t_vec(1))
     call dumpRealInterface(21,"    gyro_ti_over_te_2_in ",    gyro_ti_over_te_2_in , &
          "t_vec(2)",t_vec(2))
     call dumpRealInterface(21,"    gyro_ti_over_te_3_in ",    gyro_ti_over_te_3_in , &
          "t_vec(3)",t_vec(3))
     call dumpRealInterface(21,"    gyro_ti_over_te_4_in ",    gyro_ti_over_te_4_in , &
          "t_vec(4)",t_vec(4))
     call dumpRealInterface(21,"    gyro_ti_over_te_5_in ",    gyro_ti_over_te_5_in , &
          "t_vec(5)",t_vec(5))
     call dumpIntInterface(21,"    gyro_reintegrate_flag_in ",    gyro_reintegrate_flag_in , &
          "reintegrate_flag",reintegrate_flag)
     call dumpRealInterface(21,"    gyro_eps_dlnndr_in ",    gyro_eps_dlnndr_in , &
          "eps_dlnndr_vec(1)",eps_dlnndr_vec(1))
     call dumpRealInterface(21,"    gyro_eps_dlnndr_2_in ",    gyro_eps_dlnndr_2_in , &
          "eps_dlnndr_vec(2)",eps_dlnndr_vec(2))
     call dumpRealInterface(21,"    gyro_eps_dlnndr_3_in ",    gyro_eps_dlnndr_3_in , &
          "eps_dlnndr_vec(3)",eps_dlnndr_vec(3))
     call dumpRealInterface(21,"    gyro_eps_dlnndr_4_in ",    gyro_eps_dlnndr_4_in , &
          "eps_dlnndr_vec(4)",eps_dlnndr_vec(4))
     call dumpRealInterface(21,"    gyro_eps_dlnndr_5_in ",    gyro_eps_dlnndr_5_in , &
          "eps_dlnndr_vec(5)",eps_dlnndr_vec(5))
     call dumpRealInterface(21,"    gyro_eps_dlnndr_electron_in ",    gyro_eps_dlnndr_electron_in , &
          "eps_dlnndr_vec(0)",eps_dlnndr_vec(0))
     call dumpRealInterface(21,"    gyro_eps_dlntdr_in ",    gyro_eps_dlntdr_in , &
          "eps_dlntdr_vec(1)",eps_dlntdr_vec(1))
     call dumpRealInterface(21,"    gyro_eps_dlntdr_2_in ",    gyro_eps_dlntdr_2_in , &
          "eps_dlntdr_vec(2)",eps_dlntdr_vec(2))
     call dumpRealInterface(21,"    gyro_eps_dlntdr_3_in ",    gyro_eps_dlntdr_3_in , &
          "eps_dlntdr_vec(3)",eps_dlntdr_vec(3))
     call dumpRealInterface(21,"    gyro_eps_dlntdr_4_in ",    gyro_eps_dlntdr_4_in , &
          "eps_dlntdr_vec(4)",eps_dlntdr_vec(4))
     call dumpRealInterface(21,"    gyro_eps_dlntdr_5_in ",    gyro_eps_dlntdr_5_in , &
          "eps_dlntdr_vec(5)",eps_dlntdr_vec(5))
     call dumpRealInterface(21,"    gyro_eps_dlntdr_electron_in ",    gyro_eps_dlntdr_electron_in , &
          "eps_dlntdr_vec(0)",eps_dlntdr_vec(0))
     call dumpRealInterface(21,"    gyro_z_in ",    gyro_z_in , &
          "z_vec(1)",z_vec(1))
     call dumpRealInterface(21,"    gyro_z_2_in ",    gyro_z_2_in , &
          "z_vec(2)",z_vec(2))
     call dumpRealInterface(21,"    gyro_z_3_in ",    gyro_z_3_in , &
          "z_vec(3)",z_vec(3))
     call dumpRealInterface(21,"    gyro_z_4_in ",    gyro_z_4_in , &
          "z_vec(4)",z_vec(4))
     call dumpRealInterface(21,"    gyro_z_5_in ",    gyro_z_5_in , &
          "z_vec(5)",z_vec(5))
     call dumpRealInterface(21,"    gyro_orbit_upwind_in ",    gyro_orbit_upwind_in , &
          "orbit_upwind_vec(1)",orbit_upwind_vec(1))
     call dumpRealInterface(21,"    gyro_orbit_upwind_2_in ",    gyro_orbit_upwind_2_in , &
          "orbit_upwind_vec(2)",orbit_upwind_vec(2))
     call dumpRealInterface(21,"    gyro_orbit_upwind_3_in ",    gyro_orbit_upwind_3_in , &
          "orbit_upwind_vec(3)",orbit_upwind_vec(3))
     call dumpRealInterface(21,"    gyro_orbit_upwind_4_in ",    gyro_orbit_upwind_4_in , &
          "orbit_upwind_vec(4)",orbit_upwind_vec(4))
     call dumpRealInterface(21,"    gyro_orbit_upwind_5_in ",    gyro_orbit_upwind_5_in , &
          "orbit_upwind_vec(5)",orbit_upwind_vec(5))
     call dumpRealInterface(21,"    gyro_orbit_upwind_electron_in ",    gyro_orbit_upwind_electron_in , &
          "orbit_upwind_vec(0)",orbit_upwind_vec(0))
     call dumpRealInterface(21,"    gyro_pgamma_in ",    gyro_pgamma_in , &
          "pgamma0",pgamma0)
     call dumpRealInterface(21,"    gyro_pgamma_scale_in ",    gyro_pgamma_scale_in , &
          "pgamma0_scale",pgamma0_scale)
     call dumpRealInterface(21,"    gyro_mach_in ",    gyro_mach_in , &
          "mach0",mach0)
     call dumpRealInterface(21,"    gyro_mach_scale_in ",    gyro_mach_scale_in , &
          "mach0_scale",mach0_scale)
     call dumpIntInterface(21,"    gyro_lindiff_method_in ",    gyro_lindiff_method_in , &
          "lindiff_method",lindiff_method)
     call dumpIntInterface(21,"    gyro_trapdiff_flag_in ",    gyro_trapdiff_flag_in , &
          "trapdiff_flag",trapdiff_flag)
     call dumpIntInterface(21,"    gyro_restart_new_flag_in ",    gyro_restart_new_flag_in , &
          "restart_new_flag",restart_new_flag)
     call dumpIntInterface(21,"    gyro_restart_data_skip_in ",    gyro_restart_data_skip_in , &
          "restart_data_skip",restart_data_skip)
     call dumpIntInterface(21,"    gyro_kill_i_parallel_flag_in ",    gyro_kill_i_parallel_flag_in , &
          "kill_i_parallel_flag",kill_i_parallel_flag)
     call dumpIntInterface(21,"    gyro_kill_i_drift_flag_in ",    gyro_kill_i_drift_flag_in , &
          "kill_i_drift_flag",kill_i_drift_flag)
     call dumpIntInterface(21,"    gyro_kill_e_drift_flag_in ",    gyro_kill_e_drift_flag_in , &
          "kill_e_drift_flag",kill_e_drift_flag)
     call dumpIntInterface(21,"    gyro_kill_coll_flag_in ",    gyro_kill_coll_flag_in , &
          "kill_coll_flag",kill_coll_flag)
     call dumpRealInterface(21,"    gyro_doppler_scale_in ",    gyro_doppler_scale_in , &
          "doppler_scale",doppler_scale)
     call dumpIntInterface(21,"    gyro_nl_method_in ",    gyro_nl_method_in , &
          "nl_method",nl_method)
     call dumpIntInterface(21,"    gyro_kill_gyro_b_flag_in ",    gyro_kill_gyro_b_flag_in , &
          "kill_gyro_b_flag",kill_gyro_b_flag)
     call dumpIntInterface(21,"    gyro_velocity_output_flag_in ",    gyro_velocity_output_flag_in , &
          "velocity_output_flag",velocity_output_flag)
     call dumpIntInterface(21,"    gyro_field_r0_flag_in ",    gyro_field_r0_flag_in , &
          "field_r0_flag",field_r0_flag)
     call dumpIntInterface(21,"    gyro_field_r0_grid_in ",    gyro_field_r0_grid_in , &
          "field_r0_grid",field_r0_grid)
     call dumpRealInterface(21,"    gyro_q_scale_in ",    gyro_q_scale_in , &
          "q_scale",q_scale)
     call dumpIntInterface(21,"    gyro_dist_print_flag_in ",    gyro_dist_print_flag_in , &
          "dist_print",dist_print)
     call dumpIntInterface(21,"    gyro_nint_orb_s_in ",    gyro_nint_orb_s_in , &
          "nint_ORB_s",nint_ORB_s)
     call dumpIntInterface(21,"    gyro_nint_orb_do_in ",    gyro_nint_orb_do_in , &
          "nint_ORB_do",nint_ORB_do)
     call dumpIntInterface(21,"    gyro_udsymmetry_flag_in ",    gyro_udsymmetry_flag_in , &
          "udsymmetry_flag",udsymmetry_flag)
     call dumpIntInterface(21,"    gyro_gyro_method_in ",    gyro_gyro_method_in , &
          "gyro_method",gyro_method)
     call dumpIntInterface(21,"    gyro_sparse_method_in ",    gyro_sparse_method_in , &
          "sparse_method",sparse_method)
     call dumpIntInterface(21,"    gyro_n_mumps_max_in ",    gyro_n_mumps_max_in , &
          "n_mumps_max",n_mumps_max)
     call dumpIntInterface(21,"    gyro_n_study_in ",    gyro_n_study_in , &
          "n_study",n_study)
     call dumpRealInterface(21,"    gyro_amp_phi_study_in ",    gyro_amp_phi_study_in , &
          "amp_study",amp_study)
     call dumpRealInterface(21,"    gyro_lambda_debye_scale_in ",    gyro_lambda_debye_scale_in , &
          "lambda_debye_scale",lambda_debye_scale)
     call dumpRealInterface(21,"    gyro_lambda_debye_in ",    gyro_lambda_debye_in , &
          "lambda_debye",lambda_debye)
     call dumpIntInterface(21,"    gyro_radial_grid_offset_in ",    gyro_radial_grid_offset_in , &
          "n_x_offset",n_x_offset)
     call dumpIntInterface(21,"    gyro_theta_mult_in ",    gyro_theta_mult_in , &
          "n_theta_mult",n_theta_mult)
     call dumpIntInterface(21,"    gyro_silent_flag_in ",    gyro_silent_flag_in , &
          "silent_flag",silent_flag)
     call dumpIntInterface(21,"    gyro_nonlinear_transfer_flag_in ",    gyro_nonlinear_transfer_flag_in , &
          "nonlinear_transfer_flag",nonlinear_transfer_flag)
     call dumpRealInterface(21,"    gyro_l_x_in ",    gyro_l_x_in , &
          "l_x",l_x)
     call dumpRealInterface(21,"    gyro_l_y_in ",    gyro_l_y_in , &
          "l_y",l_y)
     call dumpIntInterface(21,"    gyro_entropy_flag_in ",    gyro_entropy_flag_in , &
          "entropy_flag",entropy_flag)
     call dumpIntInterface(21,"    gyro_ord_rbf_in ",    gyro_ord_rbf_in , &
          "ord_rbf",ord_rbf)
     call dumpIntInterface(21,"    gyro_num_equil_flag_in ",    gyro_num_equil_flag_in , &
          "num_equil_flag",num_equil_flag)
     call dumpRealInterface(21,"    gyro_zmag_in ",    gyro_zmag_in , &
          "zmag0",zmag0)
     call dumpRealInterface(21,"    gyro_dzmag_in ",    gyro_dzmag_in , &
          "dzmag0",dzmag0)
     call dumpIntInterface(21,"    gyro_output_flag_in ",    gyro_output_flag_in , &
          "output_flag",output_flag)
     call dumpRealInterface(21,"    gyro_ipccw_in ",    gyro_ipccw_in , &
          "ipccw",ipccw)
     call dumpRealInterface(21,"    gyro_btccw_in ",    gyro_btccw_in , &
          "btccw",btccw)
     call dumpIntInterface(21,"    gyro_geo_gradbcurv_flag_in ",    gyro_geo_gradbcurv_flag_in , &
          "geo_gradbcurv_flag",geo_gradbcurv_flag)
     call dumpIntInterface(21,"    gyro_geo_fastionbeta_flag_in ",    gyro_geo_fastionbeta_flag_in , &
          "geo_fastionbeta_flag",geo_fastionbeta_flag)
     call dumpRealInterface(21,"    gyro_geo_betaprime_scale_in ",    gyro_geo_betaprime_scale_in , &
          "geo_betaprime_scale",geo_betaprime_scale)
     call dumpIntInterface(21,"    gyro_poisson_z_eff_flag_in ",    gyro_poisson_z_eff_flag_in , &
          "poisson_z_eff_flag",poisson_z_eff_flag)
     call dumpIntInterface(21,"    gyro_z_eff_method_in ",    gyro_z_eff_method_in , &
          "z_eff_method",z_eff_method)
     call dumpIntInterface(21,"    gyro_gkeigen_proc_mult_in ",    gyro_gkeigen_proc_mult_in , &
          "gkeigen_proc_mult",gkeigen_proc_mult)
     call dumpIntInterface(21,"    gyro_gkeigen_method_in ",    gyro_gkeigen_method_in , &
          "gkeigen_method",gkeigen_method)
     call dumpIntInterface(21,"    gyro_gkeigen_matrixonly_in ",    gyro_gkeigen_matrixonly_in , &
          "gkeigen_matrixonly",gkeigen_matrixonly)
     call dumpIntInterface(21,"    gyro_gkeigen_mwrite_flag_in ",    gyro_gkeigen_mwrite_flag_in , &
          "gkeigen_mwrite_flag",gkeigen_mwrite_flag)
     call dumpIntInterface(21,"    gyro_gkeigen_kspace_dim_in ",    gyro_gkeigen_kspace_dim_in , &
          "gkeigen_kspace_dim",gkeigen_kspace_dim)
     call dumpIntInterface(21,"    gyro_gkeigen_n_values_in ",    gyro_gkeigen_n_values_in , &
          "gkeigen_n_values",gkeigen_n_values)
     call dumpIntInterface(21,"    gyro_gkeigen_iter_in ",    gyro_gkeigen_iter_in , &
          "gkeigen_iter",gkeigen_iter)
     call dumpRealInterface(21,"    gyro_gkeigen_tol_in ",    gyro_gkeigen_tol_in , &
          "gkeigen_tol",gkeigen_tol)
     call dumpRealInterface(21,"    gyro_gkeigen_omega_target_in ",    gyro_gkeigen_omega_target_in , &
          "gkeigen_omega_target",gkeigen_omega_target)
     call dumpRealInterface(21,"    gyro_gkeigen_gamma_target_in ",    gyro_gkeigen_gamma_target_in , &
          "gkeigen_gamma_target",gkeigen_gamma_target)
     call dumpIntInterface(21,"    gyro_linsolve_method_in ",    gyro_linsolve_method_in , &
          "linsolve_method",linsolve_method)
     call dumpIntInterface(21,"    gyro_fieldeigen_root_method_in ",    gyro_fieldeigen_root_method_in , &
          "fieldeigen_root_method",fieldeigen_root_method)
     call dumpRealInterface(21,"    gyro_fieldeigen_wr_in ",    gyro_fieldeigen_wr_in , &
          "fieldeigen_wr",fieldeigen_wr)
     call dumpRealInterface(21,"    gyro_fieldeigen_wi_in ",    gyro_fieldeigen_wi_in , &
          "fieldeigen_wi",fieldeigen_wi)
     call dumpRealInterface(21,"    gyro_fieldeigen_tol_in ",    gyro_fieldeigen_tol_in , &
          "fieldeigen_tol",fieldeigen_tol)
     call dumpIntInterface(21,"    gyro_collision_method_in ",    gyro_collision_method_in , &
          "collision_method",collision_method)
     call dumpIntInterface(21,"    gyro_io_method_in        ",    gyro_io_method_in        , &
          "io_method",io_method)
     call dumpIntInterface(21,"    gyro_time_skip_wedge_in ",    gyro_time_skip_wedge_in , &
          "time_skip_wedge",time_skip_wedge)
     call dumpIntInterface(21,"    gyro_n_torangle_wedge_in ",    gyro_n_torangle_wedge_in , &
          "n_torangle_wedge",n_torangle_wedge)
     call dumpIntInterface(21,"    gyro_n_torangle_3d_in  ",    gyro_n_torangle_3d_in  , &
          "n_torangle_3d",n_torangle_3d)
     call dumpRealInterface(21,"    gyro_theta_wedge_offset_in  ",    gyro_theta_wedge_offset_in  , &
          "theta_wedge_offset",theta_wedge_offset)
     call dumpRealInterface(21,"    gyro_theta_wedge_angle_in ",    gyro_theta_wedge_angle_in , &
          "theta_wedge_angle",theta_wedge_angle)
     call dumpRealInterface(21,"    gyro_torangle_offset_in        ",    gyro_torangle_offset_in        , &
          "torangle_offset",torangle_offset)
     close(funit)
  end if
  !
end subroutine gyro_dump_interface
!--------------------------------------------------------

 subroutine dumpIntInterface(unit1, inInterfaceVarName, inInterfaceVar, globalVarName, globalVar)

 character(*), intent(in) :: inInterfaceVarName,globalVarName  
 integer, intent(in) :: inInterfaceVar, globalVar, unit1
 character(*), parameter :: FMT1 = "(A30,A3,I5.2,A3,A30,A3,I5.2)"

 write(unit1,FMT1) inInterfaceVarName,  " = ",inInterfaceVar, " , ", globalVarName, " = ",  globalVar
 return
 end subroutine dumpIntInterface
!
!--------------------------------------------------------

 subroutine dumpRealInterface(unit1, inInterfaceVarName, inInterfaceVar, globalVarName, globalVar)

 character(*), intent(in) :: inInterfaceVarName,globalVarName  
 real, intent(in) :: inInterfaceVar, globalVar
 integer, intent(in) ::  unit1
 character(*), parameter :: FMT1 = "(A30,A3,F10.6,A3,A30,A3,F10.6)"

 write(unit1,FMT1) inInterfaceVarName,  " = ",inInterfaceVar, " , ", globalVarName, " = ",  globalVar
 return
 end subroutine dumpRealInterface
