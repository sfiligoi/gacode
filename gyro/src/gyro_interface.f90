!-------------------------------------------------------------------------
! gyro_interface.f90
!
! PURPOSE:
!  Provides interface description for GYRO.
!
! CALLING SEQUENCE:
!  call gyro_init(...)
!  set gyro_*_in variables
!  call gyro_run(...)
!  get gyro_*_out variables
!-------------------------------------------------------------------------

module gyro_interface

  implicit none

  ! Input parameters (set to default values from python/gyro_dict.py)
  integer :: gyro_radial_grid_in = 6
  integer :: gyro_orbit_grid_in = 6
  integer :: gyro_pass_grid_in = 4
  integer :: gyro_trap_grid_in = 4
  integer :: gyro_energy_grid_in = 8
  integer :: gyro_variable_egrid_flag_in = 0
  integer :: gyro_blend_grid_in = 6
  integer :: gyro_blend_fit_order_in = 3
  integer :: gyro_theta_plot_in = 1
  integer :: gyro_toroidal_min_in = 30
  integer :: gyro_toroidal_grid_in = 1
  integer :: gyro_toroidal_sep_in = 10
  integer :: gyro_toroidal_reference_in = -1
  real    :: gyro_safety_factor_in = 2.0
  real    :: gyro_shear_in = 1.0
  real    :: gyro_delta_in = 0.0
  real    :: gyro_zeta_in = 0.0
  real    :: gyro_s_delta_in = 0.0
  real    :: gyro_s_zeta_in = 0.0
  real    :: gyro_kappa_in = 1.0
  real    :: gyro_s_kappa_in = 0.0
  real    :: gyro_shift_in = 0.0
  real    :: gyro_aspect_ratio_in = 3.0
  real    :: gyro_radius_in = 0.5
  real    :: gyro_box_multiplier_in = 1.0
  real    :: gyro_rho_star_in = 0.0025
  real    :: gyro_time_step_in = 0.2
  real    :: gyro_time_max_in = 15.0
  integer :: gyro_time_skip_in = 5
  integer :: gyro_boundary_method_in = 1
  integer :: gyro_radial_derivative_band_in = 2
  integer :: gyro_radial_gyro_band_in = 3
  integer :: gyro_explicit_damp_grid_in = 8
  real    :: gyro_explicit_damp_in = 1.0
  real    :: gyro_explicit_damp_elec_in = 1.0
  integer :: gyro_debug_flag_in = 0
  integer :: gyro_nonlinear_flag_in = 0
  real    :: gyro_amp_phi_n_in = 1.0e-4
  real    :: gyro_amp_phi_0_in = 0.0
  real    :: gyro_radial_upwind_in = 1.0
  integer :: gyro_electron_method_in = 1
  integer :: gyro_radial_profile_method_in = 1
  integer :: gyro_plot_u_flag_in = 1
  integer :: gyro_plot_n_flag_in = 0
  integer :: gyro_plot_e_flag_in = 0
  integer :: gyro_plot_v_flag_in = 0
  real    :: gyro_z_eff_in = 1.0
  real    :: gyro_betae_unit_in = 0.0
  real    :: gyro_ampere_scale_in = 1.0
  integer :: gyro_n_field_in = 1
  integer :: gyro_source_flag_in = 0
  real    :: gyro_nu_source_in = 0.1
  integer :: gyro_verbose_flag_in = 0
  integer :: gyro_nonuniform_grid_flag_in = 0
  real    :: gyro_s_grid_in = 0.0
  integer :: gyro_rotation_theory_method_in = 1
  real    :: gyro_gamma_e_in = 0.0
  real    :: gyro_nu_ei_in = 0.0
  real    :: gyro_nu_ei_scale_in = 1.0
  real    :: gyro_nu_ii_scale_in = 0.0
  real    :: gyro_nu_i_krook_in = 0.0
  real    :: gyro_plot_filter_in = 1.0
  integer :: gyro_n_source_in = 1
  integer :: gyro_flat_profile_flag_in = 0
  integer :: gyro_density_method_in = 1
  integer :: gyro_integrator_method_in = 1
  real    :: gyro_energy_max_1_in = 5.0
  real    :: gyro_energy_max_2_in = 5.0
  real    :: gyro_energy_max_3_in = 5.0
  real    :: gyro_energy_max_4_in = 5.0
  real    :: gyro_energy_max_5_in = 5.0
  real    :: gyro_energy_max_electron_in = 5.0
  real    :: gyro_mu_1_in = 1.0
  real    :: gyro_mu_2_in = 1.0
  real    :: gyro_mu_3_in = 1.0
  real    :: gyro_mu_4_in = 1.0
  real    :: gyro_mu_5_in = 1.0
  real    :: gyro_mu_electron_in = 60.0
  real    :: gyro_dlnndr_in = 1.0
  real    :: gyro_dlnndr_2_in = 1.0
  real    :: gyro_dlnndr_3_in = 1.0
  real    :: gyro_dlnndr_4_in = 1.0
  real    :: gyro_dlnndr_5_in = 1.0
  real    :: gyro_dlnndr_electron_in = 1.0
  real    :: gyro_dlntdr_in = 3.0
  real    :: gyro_dlntdr_2_in = 3.0
  real    :: gyro_dlntdr_3_in = 3.0
  real    :: gyro_dlntdr_4_in = 3.0
  real    :: gyro_dlntdr_5_in = 3.0
  real    :: gyro_dlntdr_electron_in = 3.0
  real    :: gyro_ni_over_ne_in = 1.0
  real    :: gyro_ni_over_ne_2_in = 0.0
  real    :: gyro_ni_over_ne_3_in = 0.0
  real    :: gyro_ni_over_ne_4_in = 0.0
  real    :: gyro_ni_over_ne_5_in = 0.0
  real    :: gyro_ti_over_te_in = 1.0
  real    :: gyro_ti_over_te_2_in = 1.0
  real    :: gyro_ti_over_te_3_in = 1.0
  real    :: gyro_ti_over_te_4_in = 1.0
  real    :: gyro_ti_over_te_5_in = 1.0
  real    :: gyro_eps_dlnndr_in = 0.0
  real    :: gyro_eps_dlnndr_2_in = 0.0
  real    :: gyro_eps_dlnndr_3_in = 0.0
  real    :: gyro_eps_dlnndr_4_in = 0.0
  real    :: gyro_eps_dlnndr_5_in = 0.0
  real    :: gyro_eps_dlnndr_electron_in = 0.0
  real    :: gyro_eps_dlntdr_in = 0.0
  real    :: gyro_eps_dlntdr_2_in = 0.0
  real    :: gyro_eps_dlntdr_3_in = 0.0
  real    :: gyro_eps_dlntdr_4_in = 0.0
  real    :: gyro_eps_dlntdr_5_in = 0.0
  real    :: gyro_eps_dlntdr_electron_in = 0.0
  real    :: gyro_z_in = 1.0
  real    :: gyro_z_2_in = 1.0
  real    :: gyro_z_3_in = 1.0
  real    :: gyro_z_4_in = 1.0
  real    :: gyro_z_5_in = 1.0
  real    :: gyro_orbit_upwind_in = 1.0
  real    :: gyro_orbit_upwind_2_in = 1.0
  real    :: gyro_orbit_upwind_3_in = 1.0
  real    :: gyro_orbit_upwind_4_in = 1.0
  real    :: gyro_orbit_upwind_5_in = 1.0
  real    :: gyro_orbit_upwind_electron_in = 0.0
  real    :: gyro_pgamma_in = 0.0
  real    :: gyro_pgamma_scale_in = 1.0
  real    :: gyro_mach_in = 0.0
  real    :: gyro_mach_scale_in = 1.0
  integer :: gyro_lindiff_method_in = 1
  integer :: gyro_trapdiff_flag_in = 0
  integer :: gyro_restart_new_flag_in = 1
  integer :: gyro_restart_data_skip_in = 10
  integer :: gyro_kill_i_parallel_flag_in = 0
  integer :: gyro_kill_i_drift_flag_in = 0
  integer :: gyro_kill_e_drift_flag_in = 0
  integer :: gyro_kill_coll_flag_in = 1
  real    :: gyro_doppler_scale_in = 1.0
  integer :: gyro_nl_method_in = 1
  integer :: gyro_kill_gyro_b_flag_in = 0
  integer :: gyro_velocity_output_flag_in = 0
  integer :: gyro_field_r0_flag_in = 0
  integer :: gyro_field_r0_grid_in = 1
  real    :: gyro_q_scale_in = 1.0
  integer :: gyro_dist_print_flag_in = 0
  integer :: gyro_nint_orb_s_in = 64
  integer :: gyro_nint_orb_do_in = 10
  integer :: gyro_nint_geo_in = 101
  integer :: gyro_udsymmetry_flag_in = 1
  integer :: gyro_gyro_method_in = 1
  integer :: gyro_sparse_method_in = 1
  integer :: gyro_n_mumps_max_in = 1
  integer :: gyro_n_study_in = 0
  real    :: gyro_amp_phi_study_in = 0.0
  real    :: gyro_lambda_debye_scale_in = 0.0
  real    :: gyro_lambda_debye_in = 0.0
  integer :: gyro_radial_grid_offset_in = 0
  integer :: gyro_theta_mult_in = 1
  integer :: gyro_silent_flag_in = 0
  integer :: gyro_nonlinear_transfer_flag_in = 0
  real    :: gyro_l_x_in = 64.0
  real    :: gyro_l_y_in = 64.0
  integer :: gyro_entropy_flag_in = 0
  integer :: gyro_ord_rbf_in = 3
  integer :: gyro_num_equil_flag_in = 0
  real    :: gyro_zmag_in = 0.0
  real    :: gyro_dzmag_in = 0.0
  integer :: gyro_output_flag_in = 1
  real    :: gyro_ipccw_in = -1.0
  real    :: gyro_btccw_in = -1.0
  integer :: gyro_geo_gradbcurv_flag_in = 0
  real    :: gyro_geo_betaprime_scale_in = 1.0
  integer :: gyro_poisson_z_eff_flag_in = 1
  integer :: gyro_z_eff_method_in = 1
  integer :: gyro_gkeigen_proc_mult_in = 1
  integer :: gyro_gkeigen_method_in = 1
  integer :: gyro_gkeigen_matrixonly_in = 0
  integer :: gyro_gkeigen_mwrite_flag_in = 0
  integer :: gyro_gkeigen_kspace_dim_in = 300
  integer :: gyro_gkeigen_n_values_in = 10
  integer :: gyro_gkeigen_iter_in = 100
  real    :: gyro_gkeigen_tol_in = 0.0000001
  integer :: gyro_linsolve_method_in = 1
  integer :: gyro_fieldeigen_root_method_in = 1
  real    :: gyro_fieldeigen_wr_in = -0.3
  real    :: gyro_fieldeigen_wi_in = 0.2
  real    :: gyro_fieldeigen_tol_in = 1e-6
  integer :: gyro_collision_method_in = 1

  ! io related to hdf5 and diagnostics  
  integer :: gyro_iohdf5out_in = 0
  real    :: gyro_omega_exp_in = -0.0722
  real    :: gyro_zeta_offset_in = 0.0
  integer :: gyro_fine_time_skip_in = 1
  integer :: gyro_threed_time_skip_in = 1
  integer :: gyro_n_alpha_fine_in = 1
  integer :: gyro_n_alpha_threed_in = 1
  integer :: gyro_n_alpha_plot_in = 20
  real    :: gyro_theta_fine_start_in = 1
  real    :: gyro_theta_fine_end_in =1


  ! Inputs available via interface but not by INPUT
  integer :: gyro_n_fourier_geo_in = 0
  real, dimension(8,0:16) :: gyro_a_fourier_geo_in = 0.0

  ! Output parameters
  real :: gyro_elec_pflux_out = 0.0
  real :: gyro_elec_mflux_out = 0.0
  real :: gyro_elec_eflux_out = 0.0
  real :: gyro_elec_expwd_out = 0.0

  real, dimension(5) :: gyro_ion_pflux_out = 0.0
  real, dimension(5) :: gyro_ion_mflux_out = 0.0
  real, dimension(5) :: gyro_ion_eflux_out = 0.0
  real, dimension(5) :: gyro_ion_expwd_out = 0.0

  complex :: gyro_fieldeigen_omega_out
  real :: gyro_fieldeigen_error_out

contains

  ! Map GLOBAL variables to INTERFACE parameters
  subroutine map_global2interface()

    use gyro_globals

    implicit none

    gyro_radial_grid_in = n_x
    gyro_orbit_grid_in = n_theta_section
    gyro_pass_grid_in = n_pass
    gyro_trap_grid_in = n_trap
    gyro_energy_grid_in = n_energy
    gyro_variable_egrid_flag_in = variable_egrid_flag
    gyro_blend_grid_in = n_blend
    gyro_blend_fit_order_in = blend_fit_order
    gyro_theta_plot_in = n_theta_plot
    gyro_toroidal_min_in = n0
    gyro_toroidal_grid_in = n_n
    gyro_toroidal_sep_in = d_n
    gyro_toroidal_reference_in = n_ref
    gyro_safety_factor_in = q0
    gyro_shear_in = s0
    gyro_delta_in = delta0
    gyro_zeta_in = zeta0
    gyro_s_delta_in = s_delta0
    gyro_s_zeta_in = s_zeta0
    gyro_kappa_in = kappa0
    gyro_s_kappa_in = s_kappa0
    gyro_shift_in = drmaj0
    gyro_aspect_ratio_in = r_maj
    gyro_radius_in = r0
    gyro_box_multiplier_in = box_multiplier
    gyro_rho_star_in = rho_star
    gyro_time_step_in = dt
    gyro_time_max_in = time_max
    gyro_time_skip_in = time_skip
    gyro_boundary_method_in = boundary_method
    gyro_radial_derivative_band_in = m_dx
    gyro_radial_gyro_band_in = m_gyro
    gyro_explicit_damp_grid_in = n_explicit_damp
    gyro_explicit_damp_in = explicit_damp
    gyro_explicit_damp_elec_in = explicit_damp_elec
    gyro_debug_flag_in = debug_flag
    gyro_nonlinear_flag_in = nonlinear_flag
    gyro_amp_phi_n_in = amp_n
    gyro_amp_phi_0_in = amp_0
    gyro_radial_upwind_in = radial_upwind
    gyro_electron_method_in = electron_method
    gyro_radial_profile_method_in = radial_profile_method
    gyro_plot_u_flag_in = plot_u_flag
    gyro_plot_n_flag_in = plot_n_flag
    gyro_plot_e_flag_in = plot_e_flag
    gyro_plot_v_flag_in = plot_v_flag
    gyro_z_eff_in = z_eff
    gyro_betae_unit_in = betae_unit
    gyro_ampere_scale_in = ampere_scale
    gyro_n_field_in = n_field
    gyro_source_flag_in = source_flag
    gyro_nu_source_in = nu_source
    gyro_verbose_flag_in = verbose_flag
    gyro_nonuniform_grid_flag_in = nonuniform_grid_flag
    gyro_s_grid_in = s_grid
    gyro_rotation_theory_method_in = rotation_theory_method
    gyro_gamma_e_in = gamma_e
    gyro_nu_ei_in = nu_ei
    gyro_nu_ei_scale_in = nu_ei_scale
    gyro_nu_ii_scale_in = nu_ii_scale
    gyro_nu_i_krook_in = nu_i_krook
    gyro_plot_filter_in = plot_filter
    gyro_n_source_in = n_source
    gyro_flat_profile_flag_in = flat_profile_flag
    gyro_density_method_in = density_method
    gyro_integrator_method_in = integrator_method
    gyro_energy_max_1_in = energy_max_vec(1)
    gyro_energy_max_2_in = energy_max_vec(2)
    gyro_energy_max_3_in = energy_max_vec(3)
    gyro_energy_max_4_in = energy_max_vec(4)
    gyro_energy_max_5_in = energy_max_vec(5)
    gyro_energy_max_electron_in = energy_max_vec(0)
    gyro_mu_1_in = mu_vec(1)
    gyro_mu_2_in = mu_vec(2)
    gyro_mu_3_in = mu_vec(3)
    gyro_mu_4_in = mu_vec(4)
    gyro_mu_5_in = mu_vec(5)
    gyro_mu_electron_in = mu_vec(0)
    gyro_dlnndr_in = dlnndr_vec(1)
    gyro_dlnndr_2_in = dlnndr_vec(2)
    gyro_dlnndr_3_in = dlnndr_vec(3)
    gyro_dlnndr_4_in = dlnndr_vec(4)
    gyro_dlnndr_5_in = dlnndr_vec(5)
    gyro_dlnndr_electron_in = dlnndr_vec(0)
    gyro_dlntdr_in = dlntdr_vec(1)
    gyro_dlntdr_2_in = dlntdr_vec(2)
    gyro_dlntdr_3_in = dlntdr_vec(3)
    gyro_dlntdr_4_in = dlntdr_vec(4)
    gyro_dlntdr_5_in = dlntdr_vec(5)
    gyro_dlntdr_electron_in = dlntdr_vec(0)
    gyro_ni_over_ne_in = n_vec(1)
    gyro_ni_over_ne_2_in = n_vec(2)
    gyro_ni_over_ne_3_in = n_vec(3)
    gyro_ni_over_ne_4_in = n_vec(4)
    gyro_ni_over_ne_5_in = n_vec(5)
    gyro_ti_over_te_in = t_vec(1)
    gyro_ti_over_te_2_in = t_vec(2)
    gyro_ti_over_te_3_in = t_vec(3)
    gyro_ti_over_te_4_in = t_vec(4)
    gyro_ti_over_te_5_in = t_vec(5)
    gyro_eps_dlnndr_in = eps_dlnndr_vec(1)
    gyro_eps_dlnndr_2_in = eps_dlnndr_vec(2)
    gyro_eps_dlnndr_3_in = eps_dlnndr_vec(3)
    gyro_eps_dlnndr_4_in = eps_dlnndr_vec(4)
    gyro_eps_dlnndr_5_in = eps_dlnndr_vec(5)
    gyro_eps_dlnndr_electron_in = eps_dlnndr_vec(0)
    gyro_eps_dlntdr_in = eps_dlntdr_vec(1)
    gyro_eps_dlntdr_2_in = eps_dlntdr_vec(2)
    gyro_eps_dlntdr_3_in = eps_dlntdr_vec(3)
    gyro_eps_dlntdr_4_in = eps_dlntdr_vec(4)
    gyro_eps_dlntdr_5_in = eps_dlntdr_vec(5)
    gyro_eps_dlntdr_electron_in = eps_dlntdr_vec(0)
    gyro_z_in = z_vec(1)
    gyro_z_2_in = z_vec(2)
    gyro_z_3_in = z_vec(3)
    gyro_z_4_in = z_vec(4)
    gyro_z_5_in = z_vec(5)
    gyro_orbit_upwind_in = orbit_upwind_vec(1)
    gyro_orbit_upwind_2_in = orbit_upwind_vec(2)
    gyro_orbit_upwind_3_in = orbit_upwind_vec(3)
    gyro_orbit_upwind_4_in = orbit_upwind_vec(4)
    gyro_orbit_upwind_5_in = orbit_upwind_vec(5)
    gyro_orbit_upwind_electron_in = orbit_upwind_vec(0)
    gyro_pgamma_in = pgamma0
    gyro_pgamma_scale_in = pgamma0_scale
    gyro_mach_in = mach0
    gyro_mach_scale_in = mach0_scale
    gyro_lindiff_method_in = lindiff_method
    gyro_trapdiff_flag_in = trapdiff_flag
    gyro_restart_new_flag_in = restart_new_flag
    gyro_restart_data_skip_in = restart_data_skip
    gyro_kill_i_parallel_flag_in = kill_i_parallel_flag
    gyro_kill_i_drift_flag_in = kill_i_drift_flag
    gyro_kill_e_drift_flag_in = kill_e_drift_flag
    gyro_kill_coll_flag_in = kill_coll_flag
    gyro_doppler_scale_in = doppler_scale
    gyro_nl_method_in = nl_method
    gyro_kill_gyro_b_flag_in = kill_gyro_b_flag
    gyro_velocity_output_flag_in = velocity_output_flag
    gyro_field_r0_flag_in = field_r0_flag
    gyro_field_r0_grid_in = field_r0_grid
    gyro_q_scale_in = q_scale
    gyro_dist_print_flag_in = dist_print
    gyro_nint_orb_s_in = nint_ORB_s
    gyro_nint_orb_do_in = nint_ORB_do
    gyro_nint_geo_in = nint_GEO
    gyro_udsymmetry_flag_in = udsymmetry_flag
    gyro_gyro_method_in = gyro_method
    gyro_sparse_method_in = sparse_method
    gyro_n_mumps_max_in = n_mumps_max
    gyro_n_study_in = n_study
    gyro_amp_phi_study_in = amp_study
    gyro_lambda_debye_scale_in = lambda_debye_scale
    gyro_lambda_debye_in = lambda_debye
    gyro_radial_grid_offset_in = n_x_offset
    gyro_theta_mult_in = n_theta_mult
    gyro_silent_flag_in = silent_flag
    gyro_nonlinear_transfer_flag_in = nonlinear_transfer_flag
    gyro_l_x_in = l_x
    gyro_l_y_in = l_y
    gyro_entropy_flag_in = entropy_flag
    gyro_ord_rbf_in = ord_rbf
    gyro_num_equil_flag_in = num_equil_flag
    gyro_zmag_in = zmag0
    gyro_dzmag_in = dzmag0
    gyro_output_flag_in = output_flag
    gyro_ipccw_in = ipccw
    gyro_btccw_in = btccw
    gyro_geo_gradbcurv_flag_in = geo_gradbcurv_flag
    gyro_geo_betaprime_scale_in = geo_betaprime_scale
    gyro_poisson_z_eff_flag_in = poisson_z_eff_flag
    gyro_z_eff_method_in = z_eff_method
    gyro_gkeigen_proc_mult_in = gkeigen_proc_mult
    gyro_gkeigen_method_in = gkeigen_method
    gyro_gkeigen_matrixonly_in = gkeigen_matrixonly
    gyro_gkeigen_mwrite_flag_in = gkeigen_mwrite_flag
    gyro_gkeigen_kspace_dim_in = gkeigen_kspace_dim
    gyro_gkeigen_n_values_in = gkeigen_n_values
    gyro_gkeigen_iter_in = gkeigen_iter
    gyro_gkeigen_tol_in = gkeigen_tol
    gyro_linsolve_method_in = linsolve_method
    gyro_fieldeigen_root_method_in = fieldeigen_root_method
    gyro_fieldeigen_wr_in = fieldeigen_wr
    gyro_fieldeigen_wi_in = fieldeigen_wi
    gyro_fieldeigen_tol_in = fieldeigen_tol
    gyro_collision_method_in = collision_method

    gyro_iohdf5out_in        = iohdf5out
    gyro_omega_exp_in        = omega_exp
    gyro_zeta_offset_in        = zeta_offset

    gyro_fine_time_skip_in = fine_time_skip
    gyro_n_alpha_fine_in = n_alpha_fine
    gyro_n_alpha_threed_in = n_alpha_threed
    gyro_n_alpha_plot_in  = n_alpha_plot 
    gyro_theta_fine_start_in  = theta_fine_start
    gyro_theta_fine_angle_in = theta_fine_angle 


    gyro_n_fourier_geo_in = n_fourier_geo
    gyro_a_fourier_geo_in(:,:) = a_fourier_geo(:,:)

    if (debug_flag == 1 .and. i_proc == 0) then
       print *, '[map_global2interface done]'
    endif

  end subroutine map_global2interface

  ! Map INTERFACE parameters to GLOBAL variables
  subroutine map_interface2global()

    use gyro_globals

    implicit none

    n_x = gyro_radial_grid_in
    n_theta_section = gyro_orbit_grid_in
    n_pass = gyro_pass_grid_in
    n_trap = gyro_trap_grid_in
    n_energy = gyro_energy_grid_in
    variable_egrid_flag = gyro_variable_egrid_flag_in
    n_blend = gyro_blend_grid_in
    blend_fit_order = gyro_blend_fit_order_in
    n_theta_plot = gyro_theta_plot_in
    n0 = gyro_toroidal_min_in
    n_n = gyro_toroidal_grid_in
    d_n = gyro_toroidal_sep_in
    n_ref = gyro_toroidal_reference_in
    q0 = gyro_safety_factor_in
    s0 = gyro_shear_in
    delta0 = gyro_delta_in
    zeta0 = gyro_zeta_in
    s_delta0 = gyro_s_delta_in
    s_zeta0 = gyro_s_zeta_in
    kappa0 = gyro_kappa_in
    s_kappa0 = gyro_s_kappa_in
    drmaj0 = gyro_shift_in
    r_maj = gyro_aspect_ratio_in
    r0 = gyro_radius_in
    box_multiplier = gyro_box_multiplier_in
    rho_star = gyro_rho_star_in
    dt = gyro_time_step_in
    time_max = gyro_time_max_in
    time_skip = gyro_time_skip_in
    boundary_method = gyro_boundary_method_in
    m_dx = gyro_radial_derivative_band_in
    m_gyro = gyro_radial_gyro_band_in
    n_explicit_damp = gyro_explicit_damp_grid_in
    explicit_damp = gyro_explicit_damp_in
    explicit_damp_elec = gyro_explicit_damp_elec_in
    debug_flag = gyro_debug_flag_in
    nonlinear_flag = gyro_nonlinear_flag_in
    amp_n = gyro_amp_phi_n_in
    amp_0 = gyro_amp_phi_0_in
    radial_upwind = gyro_radial_upwind_in
    electron_method = gyro_electron_method_in
    radial_profile_method = gyro_radial_profile_method_in
    plot_u_flag = gyro_plot_u_flag_in
    plot_n_flag = gyro_plot_n_flag_in
    plot_e_flag = gyro_plot_e_flag_in
    plot_v_flag = gyro_plot_v_flag_in
    z_eff = gyro_z_eff_in
    betae_unit = gyro_betae_unit_in
    ampere_scale = gyro_ampere_scale_in
    n_field = gyro_n_field_in
    source_flag = gyro_source_flag_in
    nu_source = gyro_nu_source_in
    verbose_flag = gyro_verbose_flag_in
    nonuniform_grid_flag = gyro_nonuniform_grid_flag_in
    s_grid = gyro_s_grid_in
    rotation_theory_method = gyro_rotation_theory_method_in
    gamma_e = gyro_gamma_e_in
    nu_ei = gyro_nu_ei_in
    nu_ei_scale = gyro_nu_ei_scale_in
    nu_ii_scale = gyro_nu_ii_scale_in
    nu_i_krook = gyro_nu_i_krook_in
    plot_filter = gyro_plot_filter_in
    n_source = gyro_n_source_in
    flat_profile_flag = gyro_flat_profile_flag_in
    density_method = gyro_density_method_in
    integrator_method = gyro_integrator_method_in
    energy_max_vec(1) = gyro_energy_max_1_in
    energy_max_vec(2) = gyro_energy_max_2_in
    energy_max_vec(3) = gyro_energy_max_3_in
    energy_max_vec(4) = gyro_energy_max_4_in
    energy_max_vec(5) = gyro_energy_max_5_in
    energy_max_vec(0) = gyro_energy_max_electron_in
    mu_vec(1) = gyro_mu_1_in
    mu_vec(2) = gyro_mu_2_in
    mu_vec(3) = gyro_mu_3_in
    mu_vec(4) = gyro_mu_4_in
    mu_vec(5) = gyro_mu_5_in
    mu_vec(0) = gyro_mu_electron_in
    dlnndr_vec(1) = gyro_dlnndr_in
    dlnndr_vec(2) = gyro_dlnndr_2_in
    dlnndr_vec(3) = gyro_dlnndr_3_in
    dlnndr_vec(4) = gyro_dlnndr_4_in
    dlnndr_vec(5) = gyro_dlnndr_5_in
    dlnndr_vec(0) = gyro_dlnndr_electron_in
    dlntdr_vec(1) = gyro_dlntdr_in
    dlntdr_vec(2) = gyro_dlntdr_2_in
    dlntdr_vec(3) = gyro_dlntdr_3_in
    dlntdr_vec(4) = gyro_dlntdr_4_in
    dlntdr_vec(5) = gyro_dlntdr_5_in
    dlntdr_vec(0) = gyro_dlntdr_electron_in
    n_vec(1) = gyro_ni_over_ne_in
    n_vec(2) = gyro_ni_over_ne_2_in
    n_vec(3) = gyro_ni_over_ne_3_in
    n_vec(4) = gyro_ni_over_ne_4_in
    n_vec(5) = gyro_ni_over_ne_5_in
    t_vec(1) = gyro_ti_over_te_in
    t_vec(2) = gyro_ti_over_te_2_in
    t_vec(3) = gyro_ti_over_te_3_in
    t_vec(4) = gyro_ti_over_te_4_in
    t_vec(5) = gyro_ti_over_te_5_in
    eps_dlnndr_vec(1) = gyro_eps_dlnndr_in
    eps_dlnndr_vec(2) = gyro_eps_dlnndr_2_in
    eps_dlnndr_vec(3) = gyro_eps_dlnndr_3_in
    eps_dlnndr_vec(4) = gyro_eps_dlnndr_4_in
    eps_dlnndr_vec(5) = gyro_eps_dlnndr_5_in
    eps_dlnndr_vec(0) = gyro_eps_dlnndr_electron_in
    eps_dlntdr_vec(1) = gyro_eps_dlntdr_in
    eps_dlntdr_vec(2) = gyro_eps_dlntdr_2_in
    eps_dlntdr_vec(3) = gyro_eps_dlntdr_3_in
    eps_dlntdr_vec(4) = gyro_eps_dlntdr_4_in
    eps_dlntdr_vec(5) = gyro_eps_dlntdr_5_in
    eps_dlntdr_vec(0) = gyro_eps_dlntdr_electron_in
    z_vec(1) = gyro_z_in
    z_vec(2) = gyro_z_2_in
    z_vec(3) = gyro_z_3_in
    z_vec(4) = gyro_z_4_in
    z_vec(5) = gyro_z_5_in
    orbit_upwind_vec(1) = gyro_orbit_upwind_in
    orbit_upwind_vec(2) = gyro_orbit_upwind_2_in
    orbit_upwind_vec(3) = gyro_orbit_upwind_3_in
    orbit_upwind_vec(4) = gyro_orbit_upwind_4_in
    orbit_upwind_vec(5) = gyro_orbit_upwind_5_in
    orbit_upwind_vec(0) = gyro_orbit_upwind_electron_in
    pgamma0 = gyro_pgamma_in
    pgamma0_scale = gyro_pgamma_scale_in
    mach0 = gyro_mach_in
    mach0_scale = gyro_mach_scale_in
    lindiff_method = gyro_lindiff_method_in
    trapdiff_flag = gyro_trapdiff_flag_in
    restart_new_flag = gyro_restart_new_flag_in
    restart_data_skip = gyro_restart_data_skip_in
    kill_i_parallel_flag = gyro_kill_i_parallel_flag_in
    kill_i_drift_flag = gyro_kill_i_drift_flag_in
    kill_e_drift_flag = gyro_kill_e_drift_flag_in
    kill_coll_flag = gyro_kill_coll_flag_in
    doppler_scale = gyro_doppler_scale_in
    nl_method = gyro_nl_method_in
    kill_gyro_b_flag = gyro_kill_gyro_b_flag_in
    velocity_output_flag = gyro_velocity_output_flag_in
    field_r0_flag = gyro_field_r0_flag_in
    field_r0_grid = gyro_field_r0_grid_in
    q_scale = gyro_q_scale_in
    dist_print = gyro_dist_print_flag_in
    nint_ORB_s = gyro_nint_orb_s_in
    nint_ORB_do = gyro_nint_orb_do_in
    nint_GEO = gyro_nint_geo_in
    udsymmetry_flag = gyro_udsymmetry_flag_in
    gyro_method = gyro_gyro_method_in
    sparse_method = gyro_sparse_method_in
    n_mumps_max = gyro_n_mumps_max_in
    n_study = gyro_n_study_in
    amp_study = gyro_amp_phi_study_in
    lambda_debye_scale = gyro_lambda_debye_scale_in
    lambda_debye = gyro_lambda_debye_in
    n_x_offset = gyro_radial_grid_offset_in
    n_theta_mult = gyro_theta_mult_in
    silent_flag = gyro_silent_flag_in
    nonlinear_transfer_flag = gyro_nonlinear_transfer_flag_in
    l_x = gyro_l_x_in
    l_y = gyro_l_y_in
    entropy_flag = gyro_entropy_flag_in
    ord_rbf = gyro_ord_rbf_in
    num_equil_flag = gyro_num_equil_flag_in
    zmag0 = gyro_zmag_in
    dzmag0 = gyro_dzmag_in
    output_flag = gyro_output_flag_in
    ipccw = gyro_ipccw_in
    btccw = gyro_btccw_in
    geo_gradbcurv_flag = gyro_geo_gradbcurv_flag_in
    geo_betaprime_scale = gyro_geo_betaprime_scale_in
    poisson_z_eff_flag = gyro_poisson_z_eff_flag_in
    z_eff_method = gyro_z_eff_method_in
    gkeigen_proc_mult = gyro_gkeigen_proc_mult_in
    gkeigen_method = gyro_gkeigen_method_in
    gkeigen_matrixonly = gyro_gkeigen_matrixonly_in
    gkeigen_mwrite_flag = gyro_gkeigen_mwrite_flag_in
    gkeigen_kspace_dim = gyro_gkeigen_kspace_dim_in
    gkeigen_n_values = gyro_gkeigen_n_values_in
    gkeigen_iter = gyro_gkeigen_iter_in
    gkeigen_tol = gyro_gkeigen_tol_in
    linsolve_method = gyro_linsolve_method_in
    fieldeigen_root_method = gyro_fieldeigen_root_method_in
    fieldeigen_wr = gyro_fieldeigen_wr_in
    fieldeigen_wi = gyro_fieldeigen_wi_in
    fieldeigen_tol = gyro_fieldeigen_tol_in
    collision_method = gyro_collision_method_in
    
    iohdf5out      = gyro_iohdf5out_in       
    omega_exp =  gyro_omega_exp_in     
    zeta_offset = gyro_zeta_offset_in        

    fine_time_skip = gyro_fine_time_skip_in 
    threed_time_skip = gyro_threed_time_skip_in 
    n_alpha_fine = gyro_n_alpha_fine_in 
    n_alpha_threed = gyro_n_alpha_threed_in 
    n_alpha_plot = gyro_n_alpha_plot_in  
    theta_fine_start = gyro_theta_fine_start_in   
    theta_fine_end = gyro_theta_fine_end_in  


    n_fourier_geo = gyro_n_fourier_geo_in
    a_fourier_geo(:,:) = gyro_a_fourier_geo_in(:,:)

    if (debug_flag == 1 .and. i_proc == 0) then
       print *, '[map_interface2global done]'
    endif

  end subroutine map_interface2global

end module gyro_interface
