!--------------------------------------------------------------
! gyro_dump_input.f90
!
! PURPOSE:
!  Complete dump of input parameters.
!
! NOTES:
!  Necessary for debugging exactly all input variables used in a simulation.
!
! 
!---------------------------------------------------------------

 subroutine gyro_dump_input

 use gyro_globals
 implicit none
!
 integer :: funit=21
 integer :: gunit=22
 character(100) :: fname
 character(100) :: gname
!-------------------------
 if (i_proc == 0) THEN
   fname=trim(path)//"dump_inputs.txt"
   open(unit=funit,file=trim(fname),status='replace')
   gname=trim(path)//"gen_input.dat"
   open(unit=gunit,file=trim(gname),status='replace')


      call dumpIntVar(21,22,"n_x",n_x)
      call dumpIntVar(21,22,"n_theta_section",n_theta_section)
      call dumpIntVar(21,22,"n_pass",n_pass)
      call dumpIntVar(21,22,"n_trap",n_trap)
      call dumpIntVar(21,22,"n_energy",n_energy)
      call dumpIntVar(21,22,"variable_egrid_flag",variable_egrid_flag)
      call dumpIntVar(21,22,"n_blend",n_blend)
      call dumpIntVar(21,22,"blend_fit_order",blend_fit_order)
      call dumpIntVar(21,22,"n_theta_plot",n_theta_plot)
      call dumpIntVar(21,22,"n0",n0)
      call dumpIntVar(21,22,"n_n",n_n)
      call dumpIntVar(21,22,"d_n",d_n)
      call dumpIntVar(21,22,"n_ref",n_ref)
      call dumpRealVar(21,22,"q0",q0)
      call dumpRealVar(21,22,"s0",s0)
      call dumpRealVar(21,22,"delta0",delta0)
      call dumpRealVar(21,22,"zeta0",zeta0)
      call dumpRealVar(21,22,"s_delta0",s_delta0)
      call dumpRealVar(21,22,"s_zeta0",s_zeta0)
      call dumpRealVar(21,22,"kappa0",kappa0)
      call dumpRealVar(21,22,"s_kappa0",s_kappa0)
      call dumpRealVar(21,22,"drmaj0",drmaj0)
      call dumpRealVar(21,22,"r_maj",r_maj)
      call dumpRealVar(21,22,"r0",r0)
      call dumpRealVar(21,22,"box_multiplier",box_multiplier)
      call dumpRealVar(21,22,"rho_star",rho_star)
      call dumpRealVar(21,22,"dt",dt)
      call dumpRealVar(21,22,"time_max",time_max)
      call dumpIntVar(21,22,"time_skip",time_skip)
      call dumpIntVar(21,22,"boundary_method",boundary_method)
      call dumpIntVar(21,22,"m_dx",m_dx)
      call dumpIntVar(21,22,"m_gyro",m_gyro)
      call dumpIntVar(21,22,"n_explicit_damp",n_explicit_damp)
      call dumpRealVar(21,22,"explicit_damp",explicit_damp)
      call dumpRealVar(21,22,"explicit_damp_elec",explicit_damp_elec)
      call dumpIntVar(21,22,"debug_flag",debug_flag)
      call dumpIntVar(21,22,"nonlinear_flag",nonlinear_flag)
      call dumpRealVar(21,22,"amp_n",amp_n)
      call dumpRealVar(21,22,"amp_0",amp_0)
      call dumpRealVar(21,22,"radial_upwind",radial_upwind)
      call dumpIntVar(21,22,"electron_method",electron_method)
      call dumpIntVar(21,22,"radial_profile_method",radial_profile_method)
      call dumpIntVar(21,22,"plot_u_flag",plot_u_flag)
      call dumpIntVar(21,22,"plot_epar_flag",plot_epar_flag)
      call dumpIntVar(21,22,"plot_n_flag",plot_n_flag)
      call dumpIntVar(21,22,"plot_e_flag",plot_e_flag)
      call dumpIntVar(21,22,"plot_v_flag",plot_v_flag)
      call dumpRealVar(21,22,"z_eff",z_eff)
      call dumpRealVar(21,22,"betae_unit",betae_unit)
      call dumpRealVar(21,22,"ampere_scale",ampere_scale)
      call dumpIntVar(21,22,"n_field",n_field)
      call dumpIntVar(21,22,"source_flag",source_flag)
      call dumpRealVar(21,22,"nu_source",nu_source)
      call dumpIntVar(21,22,"verbose_flag",verbose_flag)
      call dumpIntVar(21,22,"nonuniform_grid_flag",nonuniform_grid_flag)
      call dumpRealVar(21,22,"s_grid",s_grid)
      call dumpIntVar(21,22,"rotation_theory_method",rotation_theory_method)
      call dumpRealVar(21,22,"gamma_e",gamma_e)
      call dumpRealVar(21,22,"nu_ei",nu_ei)
      call dumpRealVar(21,22,"nu_ei_scale",nu_ei_scale)
      call dumpRealVar(21,22,"nu_ii_scale",nu_ii_scale)
      call dumpRealVar(21,22,"nu_i_krook",nu_i_krook)
      call dumpRealVar(21,22,"plot_filter",plot_filter)
      call dumpIntVar(21,22,"n_source",n_source)
      call dumpIntVar(21,22,"flat_profile_flag",flat_profile_flag)
      call dumpIntVar(21,22,"density_method",density_method)
      call dumpIntVar(21,22,"integrator_method",integrator_method)
      call dumpRealVar(21,22," energy_max_vec(1) ", energy_max_vec(1) )
      call dumpRealVar(21,22," energy_max_vec(2) ", energy_max_vec(2) )
      call dumpRealVar(21,22," energy_max_vec(3) ", energy_max_vec(3) )
      call dumpRealVar(21,22," energy_max_vec(4) ", energy_max_vec(4) )
      call dumpRealVar(21,22," energy_max_vec(5) ", energy_max_vec(5) )
      call dumpRealVar(21,22," energy_max_vec(0) ", energy_max_vec(0) )
      call dumpRealVar(21,22," mu_vec(1) ", mu_vec(1) )
      call dumpRealVar(21,22," mu_vec(2) ", mu_vec(2) )
      call dumpRealVar(21,22," mu_vec(3) ", mu_vec(3) )
      call dumpRealVar(21,22," mu_vec(4) ", mu_vec(4) )
      call dumpRealVar(21,22," mu_vec(5) ", mu_vec(5) )
      call dumpRealVar(21,22," mu_vec(0) ", mu_vec(0) )
      call dumpRealVar(21,22," dlnndr_vec(1) ", dlnndr_vec(1) )
      call dumpRealVar(21,22," dlnndr_vec(2) ", dlnndr_vec(2) )
      call dumpRealVar(21,22," dlnndr_vec(3) ", dlnndr_vec(3) )
      call dumpRealVar(21,22," dlnndr_vec(4) ", dlnndr_vec(4) )
      call dumpRealVar(21,22," dlnndr_vec(5) ", dlnndr_vec(5) )
      call dumpRealVar(21,22," dlnndr_vec(0) ", dlnndr_vec(0) )
      call dumpRealVar(21,22," dlntdr_vec(1) ", dlntdr_vec(1) )
      call dumpRealVar(21,22," dlntdr_vec(2) ", dlntdr_vec(2) )
      call dumpRealVar(21,22," dlntdr_vec(3) ", dlntdr_vec(3) )
      call dumpRealVar(21,22," dlntdr_vec(4) ", dlntdr_vec(4) )
      call dumpRealVar(21,22," dlntdr_vec(5) ", dlntdr_vec(5) )
      call dumpRealVar(21,22," dlntdr_vec(0) ", dlntdr_vec(0) )
      call dumpRealVar(21,22," n_vec(1) ", n_vec(1) )
      call dumpRealVar(21,22," n_vec(2) ", n_vec(2) )
      call dumpRealVar(21,22," n_vec(3) ", n_vec(3) )
      call dumpRealVar(21,22," n_vec(4) ", n_vec(4) )
      call dumpRealVar(21,22," n_vec(5) ", n_vec(5) )
      call dumpRealVar(21,22," t_vec(1) ", t_vec(1) )
      call dumpRealVar(21,22," t_vec(2) ", t_vec(2) )
      call dumpRealVar(21,22," t_vec(3) ", t_vec(3) )
      call dumpRealVar(21,22," t_vec(4) ", t_vec(4) )
      call dumpRealVar(21,22," t_vec(5) ", t_vec(5) )
      call dumpRealVar(21,22," eps_dlnndr_vec(1) ", eps_dlnndr_vec(1) )
      call dumpRealVar(21,22," eps_dlnndr_vec(2) ", eps_dlnndr_vec(2) )
      call dumpRealVar(21,22," eps_dlnndr_vec(3) ", eps_dlnndr_vec(3) )
      call dumpRealVar(21,22," eps_dlnndr_vec(4) ", eps_dlnndr_vec(4) )
      call dumpRealVar(21,22," eps_dlnndr_vec(5) ", eps_dlnndr_vec(5) )
      call dumpRealVar(21,22," eps_dlnndr_vec(0) ", eps_dlnndr_vec(0) )
      call dumpRealVar(21,22," eps_dlntdr_vec(1) ", eps_dlntdr_vec(1) )
      call dumpRealVar(21,22," eps_dlntdr_vec(2) ", eps_dlntdr_vec(2) )
      call dumpRealVar(21,22," eps_dlntdr_vec(3) ", eps_dlntdr_vec(3) )
      call dumpRealVar(21,22," eps_dlntdr_vec(4) ", eps_dlntdr_vec(4) )
      call dumpRealVar(21,22," eps_dlntdr_vec(5) ", eps_dlntdr_vec(5) )
      call dumpRealVar(21,22," eps_dlntdr_vec(0) ", eps_dlntdr_vec(0) )
      call dumpRealVar(21,22," z_vec(1) ", z_vec(1) )
      call dumpRealVar(21,22," z_vec(2) ", z_vec(2) )
      call dumpRealVar(21,22," z_vec(3) ", z_vec(3) )
      call dumpRealVar(21,22," z_vec(4) ", z_vec(4) )
      call dumpRealVar(21,22," z_vec(5) ", z_vec(5) )
      call dumpRealVar(21,22," orbit_upwind_vec(1) ", orbit_upwind_vec(1) )
      call dumpRealVar(21,22," orbit_upwind_vec(2) ", orbit_upwind_vec(2) )
      call dumpRealVar(21,22," orbit_upwind_vec(3) ", orbit_upwind_vec(3) )
      call dumpRealVar(21,22," orbit_upwind_vec(4) ", orbit_upwind_vec(4) )
      call dumpRealVar(21,22," orbit_upwind_vec(5) ", orbit_upwind_vec(5) )
      call dumpRealVar(21,22," orbit_upwind_vec(0) ", orbit_upwind_vec(0) )
      call dumpRealVar(21,22,"pgamma0",pgamma0)
      call dumpRealVar(21,22,"pgamma0_scale",pgamma0_scale)
      call dumpRealVar(21,22,"mach0",mach0)
      call dumpRealVar(21,22,"mach0_scale",mach0_scale)
      call dumpIntVar(21,22,"lindiff_method",lindiff_method)
      call dumpIntVar(21,22,"trapdiff_flag",trapdiff_flag)
      call dumpIntVar(21,22,"restart_new_flag",restart_new_flag)
      call dumpIntVar(21,22,"restart_data_skip",restart_data_skip)
      call dumpIntVar(21,22,"kill_i_parallel_flag",kill_i_parallel_flag)
      call dumpIntVar(21,22,"kill_i_drift_flag",kill_i_drift_flag)
      call dumpIntVar(21,22,"kill_e_drift_flag",kill_e_drift_flag)
      call dumpIntVar(21,22,"kill_coll_flag",kill_coll_flag)
      call dumpRealVar(21,22,"doppler_scale",doppler_scale)
      call dumpIntVar(21,22,"nl_method",nl_method)
      call dumpIntVar(21,22,"kill_gyro_b_flag",kill_gyro_b_flag)
      call dumpIntVar(21,22,"velocity_output_flag",velocity_output_flag)
      call dumpIntVar(21,22,"field_r0_flag",field_r0_flag)
      call dumpIntVar(21,22,"field_r0_grid",field_r0_grid)
      call dumpRealVar(21,22,"q_scale",q_scale)
      call dumpIntVar(21,22,"dist_print",dist_print)
      call dumpIntVar(21,22,"nint_ORB_s",nint_ORB_s)
      call dumpIntVar(21,22,"nint_ORB_do",nint_ORB_do)
      call dumpIntVar(21,22,"nint_GEO",nint_GEO)
      call dumpIntVar(21,22,"udsymmetry_flag",udsymmetry_flag)
      call dumpIntVar(21,22,"gyro_method",gyro_method)
      call dumpIntVar(21,22,"sparse_method",sparse_method)
      call dumpIntVar(21,22,"n_mumps_max",n_mumps_max)
      call dumpIntVar(21,22,"n_study",n_study)
      call dumpRealVar(21,22,"amp_study",amp_study)
      call dumpRealVar(21,22,"lambda_debye_scale",lambda_debye_scale)
      call dumpRealVar(21,22,"lambda_debye",lambda_debye)
      call dumpIntVar(21,22,"n_x_offset",n_x_offset)
      call dumpIntVar(21,22,"n_theta_mult",n_theta_mult)
      call dumpIntVar(21,22,"silent_flag",silent_flag)
      call dumpIntVar(21,22,"nonlinear_transfer_flag",nonlinear_transfer_flag)
      call dumpRealVar(21,22,"l_x",l_x)
      call dumpRealVar(21,22,"l_y",l_y)
      call dumpIntVar(21,22,"entropy_flag",entropy_flag)
      call dumpIntVar(21,22,"ord_rbf",ord_rbf)
      call dumpIntVar(21,22,"num_equil_flag",num_equil_flag)
      call dumpRealVar(21,22,"zmag0",zmag0)
      call dumpRealVar(21,22,"dzmag0",dzmag0)
      call dumpIntVar(21,22,"output_flag",output_flag)
      call dumpRealVar(21,22,"ipccw",ipccw)
      call dumpRealVar(21,22,"btccw",btccw)
      call dumpIntVar(21,22,"geo_gradbcurv_flag",geo_gradbcurv_flag)
      call dumpIntVar(21,22,"geo_fastionbeta_flag",geo_fastionbeta_flag)
      call dumpRealVar(21,22,"geo_betaprime_scale",geo_betaprime_scale)
      call dumpIntVar(21,22,"poisson_z_eff_flag",poisson_z_eff_flag)
      call dumpIntVar(21,22,"z_eff_method",z_eff_method)
      call dumpIntVar(21,22,"gkeigen_proc_mult",gkeigen_proc_mult)
      call dumpIntVar(21,22,"gkeigen_method",gkeigen_method)
      call dumpIntVar(21,22,"gkeigen_matrixonly",gkeigen_matrixonly)
      call dumpIntVar(21,22,"gkeigen_mwrite_flag",gkeigen_mwrite_flag)
      call dumpIntVar(21,22,"gkeigen_kspace_dim",gkeigen_kspace_dim)
      call dumpIntVar(21,22,"gkeigen_n_values",gkeigen_n_values)
      call dumpIntVar(21,22,"gkeigen_iter",gkeigen_iter)
      call dumpRealVar(21,22,"gkeigen_tol",gkeigen_tol)
      call dumpRealVar(21,22,"gkeigen_omega_target",gkeigen_omega_target)
      call dumpRealVar(21,22,"gkeigen_gamma_target",gkeigen_gamma_target)
      call dumpIntVar(21,22,"linsolve_method",linsolve_method)
      call dumpIntVar(21,22,"fieldeigen_root_method",fieldeigen_root_method)
      call dumpRealVar(21,22,"fieldeigen_wr",fieldeigen_wr)
      call dumpRealVar(21,22,"fieldeigen_wi",fieldeigen_wi)
      call dumpRealVar(21,22,"fieldeigen_tol",fieldeigen_tol)
      call dumpIntVar(21,22,"collision_method",collision_method)
      call dumpIntVar(21,22,"io_method",io_method)
      call dumpIntVar(21,22,"time_skip_wedge",time_skip_wedge)
      call dumpIntVar(21,22,"n_torangle_wedge",n_torangle_wedge)
      call dumpIntVar(21,22,"n_torangle_3d",n_torangle_3d)
      call dumpRealVar(21,22,"theta_wedge_offset",theta_wedge_offset)
      call dumpRealVar(21,22,"theta_wedge_angle",theta_wedge_angle)
      call dumpRealVar(21,22,"torangle_offset",torangle_offset)
  close(funit)
 end if
!
 end subroutine gyro_dump_input
!--------------------------------------------------------

 subroutine dumpIntVar(unit1, unit2, inVarName, inVar)

 character(*), intent(in) :: inVarName
 integer, intent(in) :: inVar, unit1, unit2

 write(unit1,*) inVarName,  " = ",inVar
 write(unit2,*) inVar
 return
 end subroutine dumpIntVar
!
!--------------------------------------------------------

 subroutine dumpRealVar(unit1, unit2, inVarName, inVar)

 character(*), intent(in) :: inVarName
 real, intent(in) :: inVar
 integer, intent(in) :: unit1, unit2

 write(unit1,*) inVarName,  " = ",inVar
 write(unit2,*) inVar
 return
 end subroutine dumpRealVar
