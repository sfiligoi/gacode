!--------------------------------------------------------------
! gyro_input_check.f90
!
! PURPOSE:
!  Complete dump of input parameters.
!
! NOTES:
!  Necessary for debugging exactly all input variables used in a 
!  simulation.
!---------------------------------------------------------------

subroutine gyro_input_check

  use gyro_globals

  !------------------------------------
  implicit none
  !
  integer :: funit=1
  character(100) :: fname 
  !------------------------------------

  fname = trim(path)//'out.gyro.inputcheck'

  if (i_proc == 0) open(unit=funit,file=trim(fname),status='replace')

  call checkintvar(funit,"n_x",n_x)
  call checkintvar(funit,"n_theta_section",n_theta_section)
  call checkintvar(funit,"n_pass",n_pass)
  call checkintvar(funit,"n_trap",n_trap)
  call checkintvar(funit,"n_energy",n_energy)
  call checkintvar(funit,"variable_egrid_flag",variable_egrid_flag)
  call checkintvar(funit,"n_blend",n_blend)
  call checkintvar(funit,"blend_fit_order",blend_fit_order)
  call checkintvar(funit,"n_theta_plot",n_theta_plot)
  call checkintvar(funit,"n0",n0)
  call checkintvar(funit,"n_n",n_n)
  call checkintvar(funit,"d_n",d_n)
  call checkintvar(funit,"n_ref",n_ref)
  call checkrealvar(funit,"q0",q0)
  call checkrealvar(funit,"s0",s0)
  call checkrealvar(funit,"delta0",delta0)
  call checkrealvar(funit,"zeta0",zeta0)
  call checkrealvar(funit,"s_delta0",s_delta0)
  call checkrealvar(funit,"s_zeta0",s_zeta0)
  call checkrealvar(funit,"kappa0",kappa0)
  call checkrealvar(funit,"s_kappa0",s_kappa0)
  call checkrealvar(funit,"drmaj0",drmaj0)
  call checkrealvar(funit,"r_maj",r_maj)
  call checkrealvar(funit,"r0",r0)
  call checkrealvar(funit,"box_multiplier",box_multiplier)
  call checkrealvar(funit,"rho_star",rho_star)
  call checkrealvar(funit,"dt",dt)
  call checkrealvar(funit,"time_max",time_max)
  call checkintvar(funit,"time_skip",time_skip)
  call checkintvar(funit,"boundary_method",boundary_method)
  call checkintvar(funit,"m_dx",m_dx)
  call checkintvar(funit,"m_gyro",m_gyro)
  call checkintvar(funit,"n_explicit_damp",n_explicit_damp)
  call checkrealvar(funit,"explicit_damp",explicit_damp)
  call checkrealvar(funit,"explicit_damp_elec",explicit_damp_elec)
  call checkintvar(funit,"debug_flag",debug_flag)
  call checkintvar(funit,"nonlinear_flag",nonlinear_flag)
  call checkrealvar(funit,"amp_n",amp_n)
  call checkrealvar(funit,"amp_0",amp_0)
  call checkrealvar(funit,"radial_upwind",radial_upwind)
  call checkintvar(funit,"electron_method",electron_method)
  call checkintvar(funit,"radial_profile_method",radial_profile_method)
  call checkintvar(funit,"plot_u_flag",plot_u_flag)
  call checkintvar(funit,"plot_epar_flag",plot_epar_flag)
  call checkintvar(funit,"plot_n_flag",plot_n_flag)
  call checkintvar(funit,"plot_e_flag",plot_e_flag)
  call checkintvar(funit,"plot_v_flag",plot_v_flag)
  call checkrealvar(funit,"z_eff",z_eff)
  call checkrealvar(funit,"betae_unit",betae_unit)
  call checkrealvar(funit,"ampere_scale",ampere_scale)
  call checkintvar(funit,"n_field",n_field)
  call checkintvar(funit,"source_flag",source_flag)
  call checkrealvar(funit,"nu_source",nu_source)
  call checkintvar(funit,"verbose_flag",verbose_flag)
  call checkintvar(funit,"nonuniform_grid_flag",nonuniform_grid_flag)
  call checkrealvar(funit,"s_grid",s_grid)
  call checkintvar(funit,"rotation_theory_method",rotation_theory_method)
  call checkrealvar(funit,"gamma_e",gamma_e)
  call checkrealvar(funit,"nu_ei",nu_ei)
  call checkrealvar(funit,"nu_ei_scale",nu_ei_scale)
  call checkrealvar(funit,"nu_ii_scale",nu_ii_scale)
  call checkrealvar(funit,"nu_i_krook",nu_i_krook)
  call checkrealvar(funit,"plot_filter",plot_filter)
  call checkintvar(funit,"n_source",n_source)
  call checkintvar(funit,"flat_profile_flag",flat_profile_flag)
  call checkintvar(funit,"density_method",density_method)
  call checkintvar(funit,"integrator_method",integrator_method)
  call checkrealvar(funit," energy_max_vec(1) ", energy_max_vec(1) )
  call checkrealvar(funit," energy_max_vec(2) ", energy_max_vec(2) )
  call checkrealvar(funit," energy_max_vec(3) ", energy_max_vec(3) )
  call checkrealvar(funit," energy_max_vec(4) ", energy_max_vec(4) )
  call checkrealvar(funit," energy_max_vec(5) ", energy_max_vec(5) )
  call checkrealvar(funit," energy_max_vec(0) ", energy_max_vec(0) )
  call checkrealvar(funit," mu_vec(1) ", mu_vec(1) )
  call checkrealvar(funit," mu_vec(2) ", mu_vec(2) )
  call checkrealvar(funit," mu_vec(3) ", mu_vec(3) )
  call checkrealvar(funit," mu_vec(4) ", mu_vec(4) )
  call checkrealvar(funit," mu_vec(5) ", mu_vec(5) )
  call checkrealvar(funit," mu_vec(0) ", mu_vec(0) )
  call checkrealvar(funit," dlnndr_vec(1) ", dlnndr_vec(1) )
  call checkrealvar(funit," dlnndr_vec(2) ", dlnndr_vec(2) )
  call checkrealvar(funit," dlnndr_vec(3) ", dlnndr_vec(3) )
  call checkrealvar(funit," dlnndr_vec(4) ", dlnndr_vec(4) )
  call checkrealvar(funit," dlnndr_vec(5) ", dlnndr_vec(5) )
  call checkrealvar(funit," dlnndr_vec(0) ", dlnndr_vec(0) )
  call checkrealvar(funit," dlntdr_vec(1) ", dlntdr_vec(1) )
  call checkrealvar(funit," dlntdr_vec(2) ", dlntdr_vec(2) )
  call checkrealvar(funit," dlntdr_vec(3) ", dlntdr_vec(3) )
  call checkrealvar(funit," dlntdr_vec(4) ", dlntdr_vec(4) )
  call checkrealvar(funit," dlntdr_vec(5) ", dlntdr_vec(5) )
  call checkrealvar(funit," dlntdr_vec(0) ", dlntdr_vec(0) )
  call checkrealvar(funit," n_vec(1) ", n_vec(1) )
  call checkrealvar(funit," n_vec(2) ", n_vec(2) )
  call checkrealvar(funit," n_vec(3) ", n_vec(3) )
  call checkrealvar(funit," n_vec(4) ", n_vec(4) )
  call checkrealvar(funit," n_vec(5) ", n_vec(5) )
  call checkrealvar(funit," t_vec(1) ", t_vec(1) )
  call checkrealvar(funit," t_vec(2) ", t_vec(2) )
  call checkrealvar(funit," t_vec(3) ", t_vec(3) )
  call checkrealvar(funit," t_vec(4) ", t_vec(4) )
  call checkrealvar(funit," t_vec(5) ", t_vec(5) )
  call checkrealvar(funit," eps_dlnndr_vec(1) ", eps_dlnndr_vec(1) )
  call checkrealvar(funit," eps_dlnndr_vec(2) ", eps_dlnndr_vec(2) )
  call checkrealvar(funit," eps_dlnndr_vec(3) ", eps_dlnndr_vec(3) )
  call checkrealvar(funit," eps_dlnndr_vec(4) ", eps_dlnndr_vec(4) )
  call checkrealvar(funit," eps_dlnndr_vec(5) ", eps_dlnndr_vec(5) )
  call checkrealvar(funit," eps_dlnndr_vec(0) ", eps_dlnndr_vec(0) )
  call checkrealvar(funit," eps_dlntdr_vec(1) ", eps_dlntdr_vec(1) )
  call checkrealvar(funit," eps_dlntdr_vec(2) ", eps_dlntdr_vec(2) )
  call checkrealvar(funit," eps_dlntdr_vec(3) ", eps_dlntdr_vec(3) )
  call checkrealvar(funit," eps_dlntdr_vec(4) ", eps_dlntdr_vec(4) )
  call checkrealvar(funit," eps_dlntdr_vec(5) ", eps_dlntdr_vec(5) )
  call checkrealvar(funit," eps_dlntdr_vec(0) ", eps_dlntdr_vec(0) )
  call checkrealvar(funit," z_vec(1) ", z_vec(1) )
  call checkrealvar(funit," z_vec(2) ", z_vec(2) )
  call checkrealvar(funit," z_vec(3) ", z_vec(3) )
  call checkrealvar(funit," z_vec(4) ", z_vec(4) )
  call checkrealvar(funit," z_vec(5) ", z_vec(5) )
  call checkrealvar(funit," orbit_upwind_vec(1) ", orbit_upwind_vec(1) )
  call checkrealvar(funit," orbit_upwind_vec(2) ", orbit_upwind_vec(2) )
  call checkrealvar(funit," orbit_upwind_vec(3) ", orbit_upwind_vec(3) )
  call checkrealvar(funit," orbit_upwind_vec(4) ", orbit_upwind_vec(4) )
  call checkrealvar(funit," orbit_upwind_vec(5) ", orbit_upwind_vec(5) )
  call checkrealvar(funit," orbit_upwind_vec(0) ", orbit_upwind_vec(0) )
  call checkrealvar(funit,"pgamma0",pgamma0)
  call checkrealvar(funit,"pgamma0_scale",pgamma0_scale)
  call checkrealvar(funit,"mach0",mach0)
  call checkrealvar(funit,"mach0_scale",mach0_scale)
  call checkintvar(funit,"lindiff_method",lindiff_method)
  call checkintvar(funit,"trapdiff_flag",trapdiff_flag)
  call checkintvar(funit,"restart_new_flag",restart_new_flag)
  call checkintvar(funit,"restart_data_skip",restart_data_skip)
  call checkintvar(funit,"kill_i_parallel_flag",kill_i_parallel_flag)
  call checkintvar(funit,"kill_i_drift_flag",kill_i_drift_flag)
  call checkintvar(funit,"kill_e_drift_flag",kill_e_drift_flag)
  call checkintvar(funit,"kill_coll_flag",kill_coll_flag)
  call checkrealvar(funit,"doppler_scale",doppler_scale)
  call checkintvar(funit,"nl_method",nl_method)
  call checkintvar(funit,"kill_gyro_b_flag",kill_gyro_b_flag)
  call checkintvar(funit,"velocity_output_flag",velocity_output_flag)
  call checkintvar(funit,"field_r0_flag",field_r0_flag)
  call checkintvar(funit,"field_r0_grid",field_r0_grid)
  call checkrealvar(funit,"q_scale",q_scale)
  call checkintvar(funit,"dist_print",dist_print)
  call checkintvar(funit,"nint_ORB_s",nint_ORB_s)
  call checkintvar(funit,"nint_ORB_do",nint_ORB_do)
  call checkintvar(funit,"udsymmetry_flag",udsymmetry_flag)
  call checkintvar(funit,"gyro_method",gyro_method)
  call checkintvar(funit,"sparse_method",sparse_method)
  call checkintvar(funit,"n_mumps_max",n_mumps_max)
  call checkintvar(funit,"n_study",n_study)
  call checkrealvar(funit,"amp_study",amp_study)
  call checkrealvar(funit,"lambda_debye_scale",lambda_debye_scale)
  call checkrealvar(funit,"lambda_debye",lambda_debye)
  call checkintvar(funit,"n_x_offset",n_x_offset)
  call checkintvar(funit,"n_theta_mult",n_theta_mult)
  call checkintvar(funit,"silent_flag",silent_flag)
  call checkintvar(funit,"nonlinear_transfer_flag",nonlinear_transfer_flag)
  call checkrealvar(funit,"l_x",l_x)
  call checkrealvar(funit,"l_y",l_y)
  call checkintvar(funit,"entropy_flag",entropy_flag)
  call checkintvar(funit,"ord_rbf",ord_rbf)
  call checkintvar(funit,"num_equil_flag",num_equil_flag)
  call checkrealvar(funit,"zmag0",zmag0)
  call checkrealvar(funit,"dzmag0",dzmag0)
  call checkintvar(funit,"output_flag",output_flag)
  call checkrealvar(funit,"ipccw",ipccw)
  call checkrealvar(funit,"btccw",btccw)
  call checkintvar(funit,"geo_gradbcurv_flag",geo_gradbcurv_flag)
  call checkintvar(funit,"geo_fastionbeta_flag",geo_fastionbeta_flag)
  call checkrealvar(funit,"geo_betaprime_scale",geo_betaprime_scale)
  call checkintvar(funit,"poisson_z_eff_flag",poisson_z_eff_flag)
  call checkintvar(funit,"z_eff_method",z_eff_method)
  call checkintvar(funit,"gkeigen_proc_mult",gkeigen_proc_mult)
  call checkintvar(funit,"gkeigen_method",gkeigen_method)
  call checkintvar(funit,"gkeigen_matrixonly",gkeigen_matrixonly)
  call checkintvar(funit,"gkeigen_mwrite_flag",gkeigen_mwrite_flag)
  call checkintvar(funit,"gkeigen_kspace_dim",gkeigen_kspace_dim)
  call checkintvar(funit,"gkeigen_n_values",gkeigen_n_values)
  call checkintvar(funit,"gkeigen_iter",gkeigen_iter)
  call checkrealvar(funit,"gkeigen_tol",gkeigen_tol)
  call checkrealvar(funit,"gkeigen_omega_target",gkeigen_omega_target)
  call checkrealvar(funit,"gkeigen_gamma_target",gkeigen_gamma_target)
  call checkintvar(funit,"linsolve_method",linsolve_method)
  call checkintvar(funit,"fieldeigen_root_method",fieldeigen_root_method)
  call checkrealvar(funit,"fieldeigen_wr",fieldeigen_wr)
  call checkrealvar(funit,"fieldeigen_wi",fieldeigen_wi)
  call checkrealvar(funit,"fieldeigen_tol",fieldeigen_tol)
  call checkintvar(funit,"collision_method",collision_method)
  call checkintvar(funit,"io_method",io_method)
  call checkintvar(funit,"time_skip_wedge",time_skip_wedge)
  call checkintvar(funit,"n_torangle_wedge",n_torangle_wedge)
  call checkintvar(funit,"n_torangle_3d",n_torangle_3d)
  call checkrealvar(funit,"theta_wedge_offset",theta_wedge_offset)
  call checkrealvar(funit,"theta_wedge_angle",theta_wedge_angle)
  call checkrealvar(funit,"torangle_offset",torangle_offset)
  close(funit)

end subroutine gyro_input_check

subroutine checkintvar(unit,invarname,i)

  use mpi
  use gyro_globals, only : gyro_comm_world, i_err

  implicit none

  character(*), intent(in) :: invarname
  integer, intent(in) :: i,unit
  integer :: imax,imin

  call MPI_REDUCE(i,imax,1,MPI_INTEGER,MPI_MAX,0,gyro_comm_world,i_err)
  call MPI_REDUCE(i,imin,1,MPI_INTEGER,MPI_MIN,0,gyro_comm_world,i_err)

  write(unit,*) imax-imin

end subroutine checkintvar

subroutine checkrealvar(unit,invarname,x)

  use mpi
  use gyro_globals, only : gyro_comm_world, i_err

  implicit none

  character(*), intent(in) :: invarname
  integer, intent(in) :: unit
  real, intent(in) :: x
  real :: xmax,xmin

  call MPI_REDUCE(x,xmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,gyro_comm_world,i_err)
  call MPI_REDUCE(x,xmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,gyro_comm_world,i_err)

  write(unit,*) xmax-xmin

end subroutine checkrealvar


