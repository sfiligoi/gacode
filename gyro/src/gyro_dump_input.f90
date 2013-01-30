!--------------------------------------------------------------
! gyro_dump_input.f90
!
! PURPOSE:
!  Complete dump of input parameters.
!
! NOTES:
!  Necessary for debugging exactly all input variables used in a 
!  simulation.
!---------------------------------------------------------------

subroutine gyro_dump_input

  use gyro_globals

  !------------------------------------
  implicit none
  !
  integer :: funit=21
  integer :: gunit=22
  character(100) :: fname
  character(100) :: gname
  !------------------------------------

  if (i_proc == 0) then

     fname=trim(path)//'out.gyro.inputdump'
     open(unit=funit,file=trim(fname),status='replace')
     gname=trim(path)//"out.gyro.inputdumpraw"
     open(unit=gunit,file=trim(gname),status='replace')

     call dumpIntVar(funit,gunit,"n_x",n_x)
     call dumpIntVar(funit,gunit,"n_theta_section",n_theta_section)
     call dumpIntVar(funit,gunit,"n_pass",n_pass)
     call dumpIntVar(funit,gunit,"n_trap",n_trap)
     call dumpIntVar(funit,gunit,"n_energy",n_energy)
     call dumpIntVar(funit,gunit,"variable_egrid_flag",variable_egrid_flag)
     call dumpIntVar(funit,gunit,"n_blend",n_blend)
     call dumpIntVar(funit,gunit,"blend_fit_order",blend_fit_order)
     call dumpIntVar(funit,gunit,"n_theta_plot",n_theta_plot)
     call dumpIntVar(funit,gunit,"n0",n0)
     call dumpIntVar(funit,gunit,"n_n",n_n)
     call dumpIntVar(funit,gunit,"d_n",d_n)
     call dumpIntVar(funit,gunit,"n_ref",n_ref)
     call dumpRealVar(funit,gunit,"q0",q0)
     call dumpRealVar(funit,gunit,"s0",s0)
     call dumpRealVar(funit,gunit,"delta0",delta0)
     call dumpRealVar(funit,gunit,"zeta0",zeta0)
     call dumpRealVar(funit,gunit,"s_delta0",s_delta0)
     call dumpRealVar(funit,gunit,"s_zeta0",s_zeta0)
     call dumpRealVar(funit,gunit,"kappa0",kappa0)
     call dumpRealVar(funit,gunit,"s_kappa0",s_kappa0)
     call dumpRealVar(funit,gunit,"drmaj0",drmaj0)
     call dumpRealVar(funit,gunit,"r_maj",r_maj)
     call dumpRealVar(funit,gunit,"r0",r0)
     call dumpRealVar(funit,gunit,"box_multiplier",box_multiplier)
     call dumpRealVar(funit,gunit,"rho_star",rho_star)
     call dumpRealVar(funit,gunit,"dt",dt)
     call dumpRealVar(funit,gunit,"time_max",time_max)
     call dumpIntVar(funit,gunit,"time_skip",time_skip)
     call dumpIntVar(funit,gunit,"boundary_method",boundary_method)
     call dumpIntVar(funit,gunit,"m_dx",m_dx)
     call dumpIntVar(funit,gunit,"m_gyro",m_gyro)
     call dumpIntVar(funit,gunit,"n_explicit_damp",n_explicit_damp)
     call dumpRealVar(funit,gunit,"explicit_damp",explicit_damp)
     call dumpRealVar(funit,gunit,"explicit_damp_elec",explicit_damp_elec)
     call dumpIntVar(funit,gunit,"debug_flag",debug_flag)
     call dumpIntVar(funit,gunit,"nonlinear_flag",nonlinear_flag)
     call dumpRealVar(funit,gunit,"amp_n",amp_n)
     call dumpRealVar(funit,gunit,"amp_0",amp_0)
     call dumpRealVar(funit,gunit,"radial_upwind",radial_upwind)
     call dumpIntVar(funit,gunit,"electron_method",electron_method)
     call dumpIntVar(funit,gunit,"radial_profile_method",radial_profile_method)
     call dumpIntVar(funit,gunit,"plot_u_flag",plot_u_flag)
     call dumpIntVar(funit,gunit,"plot_epar_flag",plot_epar_flag)
     call dumpIntVar(funit,gunit,"plot_n_flag",plot_n_flag)
     call dumpIntVar(funit,gunit,"plot_e_flag",plot_e_flag)
     call dumpIntVar(funit,gunit,"plot_v_flag",plot_v_flag)
     call dumpRealVar(funit,gunit,"z_eff",z_eff)
     call dumpRealVar(funit,gunit,"betae_unit",betae_unit)
     call dumpRealVar(funit,gunit,"ampere_scale",ampere_scale)
     call dumpIntVar(funit,gunit,"n_field",n_field)
     call dumpIntVar(funit,gunit,"source_flag",source_flag)
     call dumpRealVar(funit,gunit,"nu_source",nu_source)
     call dumpIntVar(funit,gunit,"verbose_flag",verbose_flag)
     call dumpIntVar(funit,gunit,"nonuniform_grid_flag",nonuniform_grid_flag)
     call dumpRealVar(funit,gunit,"s_grid",s_grid)
     call dumpIntVar(funit,gunit,"rotation_theory_method",rotation_theory_method)
     call dumpRealVar(funit,gunit,"gamma_e",gamma_e)
     call dumpRealVar(funit,gunit,"nu_ei",nu_ei)
     call dumpRealVar(funit,gunit,"nu_ei_scale",nu_ei_scale)
     call dumpRealVar(funit,gunit,"nu_ii_scale",nu_ii_scale)
     call dumpRealVar(funit,gunit,"nu_i_krook",nu_i_krook)
     call dumpRealVar(funit,gunit,"plot_filter",plot_filter)
     call dumpIntVar(funit,gunit,"n_source",n_source)
     call dumpIntVar(funit,gunit,"flat_profile_flag",flat_profile_flag)
     call dumpIntVar(funit,gunit,"density_method",density_method)
     call dumpIntVar(funit,gunit,"integrator_method",integrator_method)
     call dumpRealVar(funit,gunit," energy_max_vec(1) ", energy_max_vec(1) )
     call dumpRealVar(funit,gunit," energy_max_vec(2) ", energy_max_vec(2) )
     call dumpRealVar(funit,gunit," energy_max_vec(3) ", energy_max_vec(3) )
     call dumpRealVar(funit,gunit," energy_max_vec(4) ", energy_max_vec(4) )
     call dumpRealVar(funit,gunit," energy_max_vec(5) ", energy_max_vec(5) )
     call dumpRealVar(funit,gunit," energy_max_vec(0) ", energy_max_vec(0) )
     call dumpRealVar(funit,gunit," mu_vec(1) ", mu_vec(1) )
     call dumpRealVar(funit,gunit," mu_vec(2) ", mu_vec(2) )
     call dumpRealVar(funit,gunit," mu_vec(3) ", mu_vec(3) )
     call dumpRealVar(funit,gunit," mu_vec(4) ", mu_vec(4) )
     call dumpRealVar(funit,gunit," mu_vec(5) ", mu_vec(5) )
     call dumpRealVar(funit,gunit," mu_vec(0) ", mu_vec(0) )
     call dumpRealVar(funit,gunit," dlnndr_vec(1) ", dlnndr_vec(1) )
     call dumpRealVar(funit,gunit," dlnndr_vec(2) ", dlnndr_vec(2) )
     call dumpRealVar(funit,gunit," dlnndr_vec(3) ", dlnndr_vec(3) )
     call dumpRealVar(funit,gunit," dlnndr_vec(4) ", dlnndr_vec(4) )
     call dumpRealVar(funit,gunit," dlnndr_vec(5) ", dlnndr_vec(5) )
     call dumpRealVar(funit,gunit," dlnndr_vec(0) ", dlnndr_vec(0) )
     call dumpRealVar(funit,gunit," dlntdr_vec(1) ", dlntdr_vec(1) )
     call dumpRealVar(funit,gunit," dlntdr_vec(2) ", dlntdr_vec(2) )
     call dumpRealVar(funit,gunit," dlntdr_vec(3) ", dlntdr_vec(3) )
     call dumpRealVar(funit,gunit," dlntdr_vec(4) ", dlntdr_vec(4) )
     call dumpRealVar(funit,gunit," dlntdr_vec(5) ", dlntdr_vec(5) )
     call dumpRealVar(funit,gunit," dlntdr_vec(0) ", dlntdr_vec(0) )
     call dumpRealVar(funit,gunit," n_vec(1) ", n_vec(1) )
     call dumpRealVar(funit,gunit," n_vec(2) ", n_vec(2) )
     call dumpRealVar(funit,gunit," n_vec(3) ", n_vec(3) )
     call dumpRealVar(funit,gunit," n_vec(4) ", n_vec(4) )
     call dumpRealVar(funit,gunit," n_vec(5) ", n_vec(5) )
     call dumpRealVar(funit,gunit," t_vec(1) ", t_vec(1) )
     call dumpRealVar(funit,gunit," t_vec(2) ", t_vec(2) )
     call dumpRealVar(funit,gunit," t_vec(3) ", t_vec(3) )
     call dumpRealVar(funit,gunit," t_vec(4) ", t_vec(4) )
     call dumpRealVar(funit,gunit," t_vec(5) ", t_vec(5) )
     call dumpIntVar(funit,gunit,"reintegrate_flag",reintegrate_flag)
     call dumpRealVar(funit,gunit," eps_dlnndr_vec(1) ", eps_dlnndr_vec(1) )
     call dumpRealVar(funit,gunit," eps_dlnndr_vec(2) ", eps_dlnndr_vec(2) )
     call dumpRealVar(funit,gunit," eps_dlnndr_vec(3) ", eps_dlnndr_vec(3) )
     call dumpRealVar(funit,gunit," eps_dlnndr_vec(4) ", eps_dlnndr_vec(4) )
     call dumpRealVar(funit,gunit," eps_dlnndr_vec(5) ", eps_dlnndr_vec(5) )
     call dumpRealVar(funit,gunit," eps_dlnndr_vec(0) ", eps_dlnndr_vec(0) )
     call dumpRealVar(funit,gunit," eps_dlntdr_vec(1) ", eps_dlntdr_vec(1) )
     call dumpRealVar(funit,gunit," eps_dlntdr_vec(2) ", eps_dlntdr_vec(2) )
     call dumpRealVar(funit,gunit," eps_dlntdr_vec(3) ", eps_dlntdr_vec(3) )
     call dumpRealVar(funit,gunit," eps_dlntdr_vec(4) ", eps_dlntdr_vec(4) )
     call dumpRealVar(funit,gunit," eps_dlntdr_vec(5) ", eps_dlntdr_vec(5) )
     call dumpRealVar(funit,gunit," eps_dlntdr_vec(0) ", eps_dlntdr_vec(0) )
     call dumpRealVar(funit,gunit," z_vec(1) ", z_vec(1) )
     call dumpRealVar(funit,gunit," z_vec(2) ", z_vec(2) )
     call dumpRealVar(funit,gunit," z_vec(3) ", z_vec(3) )
     call dumpRealVar(funit,gunit," z_vec(4) ", z_vec(4) )
     call dumpRealVar(funit,gunit," z_vec(5) ", z_vec(5) )
     call dumpRealVar(funit,gunit," orbit_upwind_vec(1) ", orbit_upwind_vec(1) )
     call dumpRealVar(funit,gunit," orbit_upwind_vec(2) ", orbit_upwind_vec(2) )
     call dumpRealVar(funit,gunit," orbit_upwind_vec(3) ", orbit_upwind_vec(3) )
     call dumpRealVar(funit,gunit," orbit_upwind_vec(4) ", orbit_upwind_vec(4) )
     call dumpRealVar(funit,gunit," orbit_upwind_vec(5) ", orbit_upwind_vec(5) )
     call dumpRealVar(funit,gunit," orbit_upwind_vec(0) ", orbit_upwind_vec(0) )
     call dumpRealVar(funit,gunit,"pgamma0",pgamma0)
     call dumpRealVar(funit,gunit,"pgamma0_scale",pgamma0_scale)
     call dumpRealVar(funit,gunit,"mach0",mach0)
     call dumpRealVar(funit,gunit,"mach0_scale",mach0_scale)
     call dumpIntVar(funit,gunit,"lindiff_method",lindiff_method)
     call dumpIntVar(funit,gunit,"trapdiff_flag",trapdiff_flag)
     call dumpIntVar(funit,gunit,"restart_new_flag",restart_new_flag)
     call dumpIntVar(funit,gunit,"restart_data_skip",restart_data_skip)
     call dumpIntVar(funit,gunit,"kill_i_parallel_flag",kill_i_parallel_flag)
     call dumpIntVar(funit,gunit,"kill_i_drift_flag",kill_i_drift_flag)
     call dumpIntVar(funit,gunit,"kill_e_drift_flag",kill_e_drift_flag)
     call dumpIntVar(funit,gunit,"kill_coll_flag",kill_coll_flag)
     call dumpRealVar(funit,gunit,"doppler_scale",doppler_scale)
     call dumpIntVar(funit,gunit,"nl_method",nl_method)
     call dumpIntVar(funit,gunit,"kill_gyro_b_flag",kill_gyro_b_flag)
     call dumpIntVar(funit,gunit,"velocity_output_flag",velocity_output_flag)
     call dumpRealVar(funit,gunit,"q_scale",q_scale)
     call dumpIntVar(funit,gunit,"dist_print",dist_print)
     call dumpIntVar(funit,gunit,"nint_ORB_s",nint_ORB_s)
     call dumpIntVar(funit,gunit,"nint_ORB_do",nint_ORB_do)
     call dumpIntVar(funit,gunit,"udsymmetry_flag",udsymmetry_flag)
     call dumpIntVar(funit,gunit,"gyro_method",gyro_method)
     call dumpIntVar(funit,gunit,"sparse_method",sparse_method)
     call dumpIntVar(funit,gunit,"n_mumps_max",n_mumps_max)
     call dumpIntVar(funit,gunit,"n_study",n_study)
     call dumpRealVar(funit,gunit,"amp_study",amp_study)
     call dumpRealVar(funit,gunit,"lambda_debye_scale",lambda_debye_scale)
     call dumpRealVar(funit,gunit,"lambda_debye",lambda_debye)
     call dumpIntVar(funit,gunit,"n_x_offset",n_x_offset)
     call dumpIntVar(funit,gunit,"n_theta_mult",n_theta_mult)
     call dumpIntVar(funit,gunit,"silent_flag",silent_flag)
     call dumpIntVar(funit,gunit,"nonlinear_transfer_flag",nonlinear_transfer_flag)
     call dumpRealVar(funit,gunit,"l_x",l_x)
     call dumpRealVar(funit,gunit,"l_y",l_y)
     call dumpIntVar(funit,gunit,"entropy_flag",entropy_flag)
     call dumpIntVar(funit,gunit,"ord_rbf",ord_rbf)
     call dumpIntVar(funit,gunit,"num_equil_flag",num_equil_flag)
     call dumpRealVar(funit,gunit,"zmag0",zmag0)
     call dumpRealVar(funit,gunit,"dzmag0",dzmag0)
     call dumpIntVar(funit,gunit,"output_flag",output_flag)
     call dumpRealVar(funit,gunit,"ipccw",ipccw)
     call dumpRealVar(funit,gunit,"btccw",btccw)
     call dumpIntVar(funit,gunit,"geo_gradbcurv_flag",geo_gradbcurv_flag)
     call dumpIntVar(funit,gunit,"geo_fastionbeta_flag",geo_fastionbeta_flag)
     call dumpRealVar(funit,gunit,"geo_betaprime_scale",geo_betaprime_scale)
     call dumpIntVar(funit,gunit,"poisson_z_eff_flag",poisson_z_eff_flag)
     call dumpIntVar(funit,gunit,"z_eff_method",z_eff_method)
     call dumpIntVar(funit,gunit,"gkeigen_proc_mult",gkeigen_proc_mult)
     call dumpIntVar(funit,gunit,"gkeigen_method",gkeigen_method)
     call dumpIntVar(funit,gunit,"gkeigen_matrixonly",gkeigen_matrixonly)
     call dumpIntVar(funit,gunit,"gkeigen_mwrite_flag",gkeigen_mwrite_flag)
     call dumpIntVar(funit,gunit,"gkeigen_kspace_dim",gkeigen_kspace_dim)
     call dumpIntVar(funit,gunit,"gkeigen_n_values",gkeigen_n_values)
     call dumpIntVar(funit,gunit,"gkeigen_iter",gkeigen_iter)
     call dumpRealVar(funit,gunit,"gkeigen_tol",gkeigen_tol)
     call dumpRealVar(funit,gunit,"gkeigen_omega_target",gkeigen_omega_target)
     call dumpRealVar(funit,gunit,"gkeigen_gamma_target",gkeigen_gamma_target)
     call dumpIntVar(funit,gunit,"linsolve_method",linsolve_method)
     call dumpIntVar(funit,gunit,"fieldeigen_root_method",fieldeigen_root_method)
     call dumpRealVar(funit,gunit,"fieldeigen_wr",fieldeigen_wr)
     call dumpRealVar(funit,gunit,"fieldeigen_wi",fieldeigen_wi)
     call dumpRealVar(funit,gunit,"fieldeigen_tol",fieldeigen_tol)
     call dumpIntVar(funit,gunit,"io_method",io_method)
     call dumpIntVar(funit,gunit,"time_skip_wedge",time_skip_wedge)
     call dumpIntVar(funit,gunit,"n_torangle_wedge",n_torangle_wedge)
     call dumpIntVar(funit,gunit,"n_torangle_3d",n_torangle_3d)
     call dumpRealVar(funit,gunit,"theta_wedge_offset",theta_wedge_offset)
     call dumpRealVar(funit,gunit,"theta_wedge_angle",theta_wedge_angle)
     call dumpRealVar(funit,gunit,"torangle_offset",torangle_offset)
     close(funit)
     close(gunit)
  endif
 
end subroutine gyro_dump_input


subroutine dumpIntVar(unit1, unit2, inVarName, inVar)

  character(*), intent(in) :: inVarName
  integer, intent(in) :: inVar, unit1, unit2

  write(unit1,*) inVarName," = ",inVar
  !write(unit2,*) inVar

end subroutine dumpIntVar


subroutine dumpRealVar(unit1, unit2, inVarName, inVar)

  character(*), intent(in) :: inVarName
  real, intent(in) :: inVar
  integer, intent(in) :: unit1, unit2

  write(unit1,*) inVarName," = ",inVar
  !write(unit2,*) inVar

end subroutine dumpRealVar
