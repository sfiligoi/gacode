!-----------------------------------------------------
! gyro_globals.f90
!
! PURPOSE:
!  Fundamental module containing most of the shared 
!  GYRO variables and arrays
!-----------------------------------------------------

module gyro_globals

  integer :: uflag

  integer :: gyrotest_flag
  integer :: lskipinit

  integer :: gyro_exit_status
  character(len=80) :: gyro_exit_message

  !---------------------------------------------------------
  ! Set this flag to unity if subgrouping
  ! is to be used (suggested).
  !
  integer, parameter :: USE_SUBGROUP = 1
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Unit for restart files:
  !
  integer, parameter :: io_restart = 2
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Length of entropy vector, entropy(:)
  !
  integer, parameter :: n_entro=5
  !------------------------------------------------------

  !------------------------------------------------------
  ! Run file:
  !
  character(len=80) :: runfile = 'run.out'
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Newline characters:
  !
  character(len=*), parameter :: new_line = char(13)//char(10)
  character(len=*), parameter :: separator = &
       '-----------------------------------------------'
  !---------------------------------------------------------

  !--------------------------------------------------
  ! Proton mass in kg:
  !
  real, parameter :: kg_proton = 1.6726e-27
  !--------------------------------------------------

  !---------------------------------------------------------
  ! Restart files:
  !
  ! in read_restart:
  !
  character(len=9), dimension(0:1), parameter :: &
       file_restart = (/ &
       "RESTART_0",&
       "RESTART_1" &
       /)
  !
  ! tag files for read_restart:
  !
  character(len=13), dimension(0:1), parameter :: &
       file_tag_restart = (/ &
       "RESTART_tag_0",&
       "RESTART_tag_1" &
       /)
  !
  ! in read_eq_restart:
  !
  character(len=12), dimension(0:1), parameter :: &
       file_eq_restart = (/ &
       "RESTART_eq_0",&
       "RESTART_eq_1" &
       /)
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Algorithm selectors:
  !
  ! (a) methods (values 1,2,...)
  !
  integer :: time_method
  integer :: boundary_method 
  integer :: electron_method
  integer :: radial_profile_method
  integer :: geometry_method
  integer :: nonuniform_grid_flag
  integer :: nonlinear_transfer_flag
  integer :: density_method
  integer :: integrator_method
  integer :: nl_method
  integer :: rotation_theory_method
  integer :: lindiff_method
  integer :: gyro_method
  integer :: sparse_method
  integer :: transport_method
  integer :: linsolve_method 
  integer :: collision_method
  integer :: fieldeigen_root_method
  !
  ! (b) flags (0 or 1)
  !
  integer :: nonlinear_flag 
  integer :: collision_flag
  integer :: krook_flag
  integer :: verbose_flag
  integer :: debug_flag
  integer :: flat_profile_flag
  integer :: trapdiff_flag
  integer :: kill_i_parallel_flag
  integer :: kill_i_drift_flag
  integer :: kill_e_drift_flag
  integer :: kill_coll_flag
  integer :: kill_gyro_b_flag
  integer :: velocity_output_flag
  integer :: geo_array_print_flag
  integer :: source_flag
  integer :: field_r0_flag
  integer :: dist_print
  integer :: udsymmetry_flag
  integer :: silent_flag
  integer :: eparallel_plot_flag
  integer :: entropy_flag
  integer :: num_equil_flag
  integer :: gkeigen_proc_mult
  integer :: gkeigen_method
  integer :: gkeigen_matrixonly
  integer :: gkeigen_mwrite_flag
  integer :: plot_u_flag
  integer :: plot_n_flag
  integer :: plot_e_flag
  integer :: plot_v_flag
  integer :: poisson_z_eff_flag
  integer :: z_eff_method
  integer :: variable_egrid_flag
  integer :: geo_gradbcurv_flag
  integer :: fakefield_flag
  !---------------------------------------------------------

  !-----------------------------------------------------------------------------------
  ! Gyrokinetic eigensolver (GKEIGEN) parameters: 
  !
  ! Number of eigenvalues to find
  integer :: gkeigen_n_values
  !
  ! Maximum Krylov subspace dimension
  integer :: gkeigen_kspace_dim      
  !
  ! Maximum number of iterations of Rayleigh-Ritz algorithm
  integer :: gkeigen_iter   
  !
  ! Total length of state vectors
  integer :: h_length     
  !
  ! Local length and width of state vector
  integer :: h_length_loc
  integer :: h_width_loc 
  integer :: h_length_block
  integer :: h_length_block_t
  integer :: seq_length  ! Number of elements in passed matrix block
  integer :: seq_length_t
  !
  ! Absolute tolerance in eigenvector solution
  real :: gkeigen_tol       
  !
  ! For added parallelization layer (gkeigen_proc_mult > 1)
  integer :: gkeigen_j_set
  integer :: j_proc_tot
  integer :: n_proc_tot
  integer :: gkeigen_transpose_flag
  !
  !-----------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------
  ! Field eigensolver (fieldeigen) parameters
  !
  real :: fieldeigen_wr
  real :: fieldeigen_wi
  real :: fieldeigen_tol
  !-----------------------------------------------------------------------------------

  !---------------------------------------------------------
  ! Grid dimensions:
  !
  integer :: n_x
  integer :: n_x_offset
  integer :: n_theta_section
  integer :: n_blend
  integer :: n_interval
  integer :: n_theta_plot
  integer :: n_theta_mult
  integer :: n_theta_int
  integer :: blend_fit_order
  integer :: n_stack 
  integer :: n_field
  integer :: nint_ORB_s
  integer :: nint_ORB_do
  integer :: nint_GEO
  integer :: n_mumps_max
  integer :: n_study
  !
  real, dimension(:), allocatable :: energy_max
  real, dimension(0:10) :: energy_max_vec
  real :: box_multiplier
  real :: l_x
  real :: l_y
  !
  integer :: n_ref
  integer :: n_lambda
  integer :: n_pass
  integer :: n_trap
  integer :: n_energy
  integer :: n_n
  integer :: n_n_1
  integer :: n0
  integer :: d_n
  integer :: n_max
  integer :: n_moment
  integer :: p_moment
  integer :: n_substep
  !
  integer :: n_lump
  integer :: field_r0_grid
  integer :: ord_rbf
  integer :: n_coll
  integer :: n_rbf
  !
  ! Species 
  ! - we always have n_spec = n_ion+1
  ! - with adiabatic electrons, n_kinetic = n_ion
  ! - with kinetic electrons, n_kinetic = n_spec
  ! - indx_e is the electron index
  !
  integer :: n_spec
  integer :: n_ion
  integer :: n_kinetic
  integer :: n_gk
  integer :: indx_e
  integer :: n_grid_exp
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Loop indices (we try to be consistent)
  !
  ! Radius
  integer :: i
  integer :: ip
  integer :: i_diff
  !
  ! Energy
  integer :: ie
  integer :: k 
  !
  ! Pitch angle
  integer :: ck
  !
  ! Blending index
  integer :: j
  integer :: jp
  !
  ! Orbit time
  integer :: m
  integer :: m0
  !
  ! Interpolated theta index
  integer :: j_int
  !
  ! Species
  integer :: is
  !
  ! Toroidal mode number
  integer :: in
  integer :: in_1
  !
  ! Field
  integer :: ix
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Field matrices/dimensions for maxwell solve:
  !
  integer :: n_maxwell
  integer :: n_maxwell_row
  !
  complex, allocatable, target, dimension(:) :: m_maxwell
  integer, allocatable, target, dimension(:) :: indx_maxwell
  !
  ! Field matrices/dimensions for collision-Ampere solve:
  !
  integer :: n_ampere
  integer :: n_ampere_row
  !
  complex, allocatable, target, dimension(:) :: m_ampere
  integer, allocatable, target, dimension(:) :: indx_ampere
  !
  integer :: n_poisson
  integer :: n_poisson_row
  !
  complex, allocatable, target, dimension(:) :: m_poisson
  integer, allocatable, target, dimension(:) :: indx_poisson
  !
  integer :: n_poissonaperp
  integer :: n_poissonaperp_row
  !
  complex, allocatable, target, dimension(:) :: m_poissonaperp
  integer, allocatable, target, dimension(:) :: indx_poissonaperp
  !
  complex, dimension(:,:,:), allocatable :: coll_vel
  complex, dimension(:,:,:), allocatable :: coll_vel_perp1, coll_vel_perp2
  !
  complex, pointer, dimension(:) :: m_sparse
  integer, pointer, dimension(:) :: indx_sparse
  !
  integer, dimension(4) :: lvalue
  integer, dimension(4) :: lindx
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Maxwell blending-function expansion coefficients:
  !
  complex, allocatable, dimension(:,:,:) :: field_blend
  complex, allocatable, dimension(:,:,:) :: field_blend_old
  complex, allocatable, dimension(:,:,:) :: field_blend_old2
  complex, allocatable, dimension(:,:,:) :: field_blend_dot
  !--------------------------------------------------------

  !---------------------------------------------------------
  ! Maxwell blending-function matrices:
  !
  integer, allocatable, dimension(:,:,:) :: ff_mm_piv
  complex, allocatable, dimension(:,:,:,:) :: ff_mm
  integer, allocatable, dimension(:,:) :: ff2_mm_piv
  complex, allocatable, dimension(:,:,:) :: ff2_mm
  !
  complex, allocatable, dimension(:,:,:,:) :: ap_mm
  complex, allocatable, dimension(:,:,:,:) :: aa_mm
  complex, allocatable, dimension(:,:,:,:) :: ab_mm
  complex, allocatable, dimension(:,:,:,:) :: abp_mm

  !
  complex, dimension(:,:,:,:), allocatable :: cs_blend
  complex, dimension(:,:,:,:), allocatable :: c_blend
  !
  complex, dimension(:,:,:,:), allocatable :: cs_blend_prime
  !
  complex, dimension(:,:,:), allocatable :: blend_plot
  complex, dimension(:,:,:), allocatable :: blend_prime_plot
  complex, dimension(:,:), allocatable :: blend_r0_plot
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Timestep parameters:
  !
  integer :: step
  integer :: nstep
  integer :: data_step
  integer :: time_skip
  integer :: alltime_index
  !
  real :: time_max
  real :: freq_tol
  real :: freq_err
  !
  real :: plot_filter
  real :: dt
  real :: t_current
  !
  real, dimension(:), allocatable :: w_time
  !
  complex, dimension(2) :: freq_n
  ! 
  real, dimension(:), allocatable :: time_error
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Restart parameters:
  !
  ! (input)
  !
  integer :: restart_method
  integer :: restart_new_flag
  integer :: restart_data_skip
  integer :: eigensolve_restart_flag
  !
  ! (working)
  !
  integer :: i_restart
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Timers:
  !
  integer :: clock_count
  integer :: clock_rate
  integer :: clock_max
  !
  real :: elapsed_time
  !
  real :: CPU_0
  real :: CPU_1
  real :: CPU_2
  real :: CPU_2a
  real :: CPU_3
  real :: CPU_4
  real :: CPU_5
  real :: CPU_6
  real :: CPU_7
  real :: CPU_C_in
  real :: CPU_C_out
  real :: CPU_C
  real :: CPU_Ct_in
  real :: CPU_Ct_out
  real :: CPU_Ct
  real :: CPU_NL_in
  real :: CPU_NL_out
  real :: CPU_NL
  real :: CPU_NLt_in
  real :: CPU_NLt_out
  real :: CPU_NLt
  real :: CPU_RHS_in
  real :: CPU_RHS_out
  real :: CPU_RHS
  real :: CPU_diag_in
  real :: CPU_diag_mid
  real :: CPU_diag_outp
  real :: CPU_diag_out
  real :: CPU_diag_a
  real :: CPU_diag_b
  real :: CPU_field_in
  real :: CPU_field_out
  real :: CPU_field
  real :: CPU_field2_in
  real :: CPU_field2_out
  real :: CPU_field2
  real :: CPU_interp
  real :: CPU_interp_in
  real :: CPU_interp_out
  real :: CPU_ts
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Radial boundary parameters:
  !
  integer, dimension(:), allocatable :: i_cyc
  integer, dimension(:), allocatable :: i_loop
  !
  integer :: i1_buffer
  integer :: i2_buffer
  !
  integer :: n_explicit_damp
  !
  real :: explicit_damp
  real :: explicit_damp_elec
  real, dimension(:,:), allocatable :: explicit_damp_vec
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Radial grid parameters
  !
  integer :: m_dx
  integer :: mg_dx
  integer :: m_gyro
  integer :: i_dx
  integer :: ig_dx
  integer :: i_gyro
  integer :: ir_norm
  !
  integer :: n_source
  !
  real :: d_x
  real :: radial_upwind
  !
  complex, dimension(:,:), allocatable :: cr
  complex, dimension(:,:), allocatable :: cri
  !
  ! Source basis functions, coupling matrix and pivots:
  !
  real, dimension(:,:), allocatable :: b_src
  real, dimension(:,:), allocatable :: m_src
  integer, dimension(:), allocatable :: src_piv
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Radial operators:
  !
  complex, dimension(:), allocatable :: w_g0
  complex, dimension(:), allocatable :: w_d0
  complex, dimension(:), allocatable :: w_gd0
  complex, dimension(:), allocatable :: w_d1
  complex, dimension(:), allocatable :: w_d2
  !
  real, dimension(:), allocatable :: s_d1
  real, dimension(:,:), allocatable :: gyro_trace
  !
  complex, dimension(:,:,:,:,:), allocatable :: w_gyro
  complex, dimension(:,:,:,:,:), allocatable :: w_gyro_rot
  complex, dimension(:,:,:,:,:), allocatable :: w_gyro_aperp
  !
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Toroidal grid parameters:
  !
  integer, dimension(:), allocatable :: n
  integer, dimension(:), allocatable :: n_1 
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Orbit grid arrays and parameters:
  !
  integer, parameter :: n_class = 2
  !
  integer, dimension(n_class) :: n_theta
  integer, dimension(n_class) :: n_tau
  !
  real, dimension(n_class)    :: d_tau
  real, dimension(n_class,2)  :: sigma
  real, dimension(n_class,2)  :: sigma_tau
  !
  integer, dimension(:), allocatable :: class
  !
  real, dimension(:,:,:), allocatable :: theta_t
  real, dimension(:,:,:), allocatable :: b0_t
  real, dimension(:,:,:), allocatable :: g_theta_t
  real, dimension(:,:,:), allocatable :: grad_r_t
  real, dimension(:,:,:), allocatable :: qrat_t
  real, dimension(:,:,:), allocatable :: cos_t
  real, dimension(:,:,:), allocatable :: cos_p_t
  real, dimension(:,:,:), allocatable :: captheta_t
  real, dimension(:,:,:), allocatable :: sin_t
  real, dimension(:,:,:), allocatable :: usin_t
  real, dimension(:,:,:), allocatable :: ucos_t
  !
  real, dimension(:,:,:), allocatable :: bt_t
  real, dimension(:,:,:), allocatable :: bp_t
  real, dimension(:,:,:), allocatable :: bigr_t
  !
  real, dimension(:,:,:), allocatable :: tau
  !
  integer, dimension(:,:,:), allocatable :: m_cyc
  complex, dimension(:,:,:,:), allocatable :: p_cyc
  !
  integer, dimension(:,:), allocatable :: m_phys
  integer, dimension(:,:), allocatable :: p_phys
  !
  real, dimension(:,:), allocatable :: omega
  !
  integer, dimension(:,:,:), allocatable :: m_map
  !
  real, dimension(:,:), allocatable :: c_fluxave
  !
  real :: d_theta
  !
  ! Implicit advection arrays and parameters
  !
  real :: a_SDIRK
  !
  complex, dimension(:,:,:,:), allocatable :: o_f
  complex, dimension(:,:,:,:), allocatable :: o_fv
  complex, dimension(:,:,:,:), allocatable :: imp
  complex, dimension(:,:,:,:), allocatable :: o_advect
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Lambda/energy grid variables:
  !
  real, allocatable, dimension(:,:) :: energy
  real, allocatable, dimension(:,:) :: w_energy
  !
  real, allocatable, dimension(:,:) :: lambda
  real, allocatable, dimension(:) :: s_lambda
  real, allocatable, dimension(:,:) :: w_lambda
  !
  real, allocatable, dimension(:,:,:,:) :: w_p
  !
  real, allocatable, dimension(:) :: lambda_tp
  real, allocatable, dimension(:) :: lambda_max
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Equilibrium specification parameters:
  !
  real :: a_meters
  !
  real :: q0
  real :: s0
  !
  real :: ipccw
  real :: btccw
  !
  real :: r_norm
  real :: q_norm
  real :: shat_norm
  real :: rhos_norm
  real :: csda_norm
  real :: tem_norm
  real :: den_norm
  real :: r_maj_norm
  real :: b_unit_norm
  real :: betae_unit_norm
  real :: betai_unit_norm
  !   
  ! Local profile parameters
  ! 
  real :: kappa0
  real :: s_kappa0
  real :: delta0
  real :: zeta0
  real :: s_delta0
  real :: s_zeta0
  real :: drmaj0
  real :: zmag0
  real :: dzmag0
  real :: alpha_mhd0
  real :: pgamma0
  real :: mach0
  real :: r_maj
  real :: r0
  real :: x_length
  real :: rho_star
  real :: z_eff
  real :: betae_unit
  real :: ampere_scale
  real :: geo_betaprime_scale
  real :: gamma_e
  real :: lambda_debye_scale
  real :: lambda_debye
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Local profile input parameters
  !  
  ! ln_vec -> L_n
  ! lt_vec -> L_T
  !  t_vec -> T
  !  n_vec -> n
  !  z_vec -> Z
  !
  ! index 0 -> electron
  ! index 1 -> ion 1
  ! etc.
  !
  real :: dlnndr_vec(0:10)
  real :: dlntdr_vec(0:10)
  real :: eps_dlnndr_vec(0:10)
  real :: eps_dlntdr_vec(0:10)
  real :: n_vec(0:10)
  real :: t_vec(0:10)
  real :: z_vec(0:10)
  real :: orbit_upwind_vec(0:10) 
  !
  ! mu_i = sqrt(m(1)/m(i))
  !
  real ::  mu_vec(0:10)
  !
  real, dimension(:,:), allocatable :: krho_i
  real, dimension(:), allocatable :: z
  real, dimension(:), allocatable :: mu
  !
  real :: nu_ei
  real :: nu_ei_scale
  real :: nu_ii_scale
  real :: nu_i_krook
  !
  real :: s_grid
  !
  ! Scales experimental ExB
  real :: doppler_scale
  ! 
  ! Scales experimental q-profile but preserves s
  real :: q_scale 
  !
  ! Scales experimental pgamma  
  real :: pgamma0_scale
  !
  ! Scales experimental  mach 
  real :: mach0_scale
  !
  ! Geometry Fourier coefficients
  integer :: n_fourier_geo
  real, dimension(8,0:16) :: a_fourier_geo
  !---------------------------------------------------------

  !---------------------------------------------
  ! Profile functions:  
  !
  ! These are "special" values which retain spatial 
  ! variation for a flux-tube run:
  !
  real, dimension(:), allocatable :: r 
  real, dimension(:), allocatable :: q
  !
  ! These get flattened like all other _s profiles
  ! in a flux-tube run
  !
  real, dimension(:), allocatable :: q_s 
  real, dimension(:), allocatable :: r_s 
  !
  real, dimension(:), allocatable :: rmaj_s
  real, dimension(:), allocatable :: rhogrid_s 
  real, dimension(:), allocatable :: kappa_s 
  real, dimension(:), allocatable :: s_kappa_s 
  real, dimension(:), allocatable :: delta_s  
  real, dimension(:), allocatable :: zeta_s  
  real, dimension(:), allocatable :: s_delta_s 
  real, dimension(:), allocatable :: s_zeta_s 
  real, dimension(:), allocatable :: z_eff_s 
  real, dimension(:), allocatable :: b_unit_s
  real, dimension(:), allocatable :: rhosda_s 
  real, dimension(:), allocatable :: csda_s  
  real, dimension(:), allocatable :: drmaj_s 
  real, dimension(:), allocatable :: zmag_s
  real, dimension(:), allocatable :: dzmag_s 
  real, dimension(:), allocatable :: shat_s 
  real, dimension(:), allocatable :: beta_unit_s 
  real, dimension(:), allocatable :: w0_s
  real, dimension(:), allocatable :: w0p_s
  real, dimension(:), allocatable :: gamma_e_s
  real, dimension(:), allocatable :: gamma_p_s
  real, dimension(:), allocatable :: mach_s
  !
  real, dimension(:), allocatable :: omega_eb_s
  real, dimension(:), allocatable :: dlnpdr_s
  real, dimension(:), allocatable :: beta_star_s
  !
  real, dimension(:,:), allocatable :: den_s 
  real, dimension(:,:), allocatable :: tem_s
  real, dimension(:,:), allocatable :: dlnndr_s 
  real, dimension(:,:), allocatable :: dlntdr_s 
  real, dimension(:,:), allocatable :: alpha_s
  real, dimension(:,:), allocatable :: nu_s   
  real, dimension(:,:), allocatable :: pr_s
  !
  complex, dimension(:,:), allocatable :: phase
  real, dimension(:), allocatable :: angp
  !
  ! nonuniform_grid_flag variables:
  !
  real, dimension(:), allocatable :: r_e
  real, dimension(:), allocatable :: dr_eodr
  !
  ! General geometry Fourier coefficients
  !
  real, dimension(:,:,:), allocatable :: a_fourier_geo_s
  !----------------------------------------------

  !---------------------------------------------------------
  ! Field and field coefficients:
  !
  real, dimension(:), allocatable :: phi_squared
  real, dimension(:,:), allocatable :: g_squared
  real, dimension(:), allocatable :: phi_fluxave
  real, dimension(:), allocatable :: a_fluxave
  real, dimension(:), allocatable :: aperp_fluxave
  !
  complex, dimension(:,:,:), allocatable :: phi
  complex, dimension(:,:,:), allocatable :: phi_plot
  complex, dimension(:,:), allocatable :: field_r0_plot
  real, dimension(:,:), allocatable :: ave_phi
  !
  complex, dimension(:,:,:,:), allocatable :: moments_plot
  real, dimension(:,:,:), allocatable :: moments_zero_plot
  !
  real, dimension(:,:), allocatable :: b0_plot
  real, dimension(:,:), allocatable :: g_theta_plot   
  !
  ! Work arrays
  !
  complex, allocatable, dimension(:,:) :: vel_sum_p
  complex, allocatable, dimension(:,:) :: vel_sum_a
  complex, allocatable, dimension(:,:) :: vel_sum_aperp
  ! 
  ! Normalization
  !
  complex :: balloon_renorm
  !
  ! Initial amplitudes (see make_initial_h)
  !
  real :: amp_n
  real :: amp_0
  real :: amp_study
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Gyroaverage of h_i
  !
  complex, dimension(:,:,:,:), allocatable :: gyro_h
  complex, dimension(:,:,:,:), allocatable :: gyro_h_aperp
  !
  ! Phi and A mapped onto the orbit grid (i.e., as
  ! functions of tau).
  !
  complex, dimension(:,:,:,:), allocatable :: field_tau
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Arrays connected with distribution function:
  !
  complex, dimension(:,:,:,:), allocatable :: h_err
  complex, dimension(:,:,:,:), allocatable :: h_0
  complex, dimension(:,:,:,:), allocatable :: h_old
  complex, dimension(:,:,:,:), allocatable :: h
  complex, dimension(:,:,:,:), allocatable :: h_cap
  complex, dimension(:,:,:,:), allocatable :: h_cap_old
  complex, dimension(:,:,:,:), allocatable :: h_cap_old2
  complex, dimension(:,:,:,:), allocatable :: h_cap_dot
  complex, dimension(:,:,:,:), allocatable :: RHS

  complex, dimension(:,:,:,:), allocatable :: f_store
  complex, dimension(:,:,:,:), allocatable :: p_store
  !
  complex, dimension(:,:,:,:), allocatable :: h_M
  complex, dimension(:,:,:), allocatable :: h_C
  complex, dimension(:,:,:), allocatable :: f_coll
  complex, dimension(:,:,:), allocatable :: fb_coll
  !
  ! Various averages of field quantities:
  !
  complex, dimension(:,:,:,:,:), allocatable :: gyro_uv
  complex, dimension(:,:,:,:,:), allocatable :: kyro_uv
  complex, dimension(:,:,:,:,:), allocatable :: gyro_uv_old2
  complex, dimension(:,:,:,:,:), allocatable :: gyro_uv_old
  complex, dimension(:,:,:,:,:), allocatable :: gyro_uv_dot
  complex, dimension(:,:,:,:), allocatable :: gyro_u
  complex, dimension(:,:,:,:), allocatable :: gyro_u_M
  !
  complex, dimension(:,:,:), allocatable :: RHS_krook
  !
  ! Collision operator
  !
  real :: condition_number
  real, dimension(:,:,:,:), allocatable :: d_rbf
  real, dimension(:,:), allocatable :: d1_rbf
  complex, dimension(:,:,:,:), allocatable :: nu_op
  !
  ! For ebelli collision method (generally not used)
  !
  complex, dimension(:,:,:,:), allocatable :: h_C_all
  real, dimension(:,:,:,:), allocatable   ::  d_rbf_lorentz
  real, dimension(:,:,:,:,:), allocatable ::  d_rbf_rs
  real, dimension(:,:,:,:), allocatable   ::  d_rbf_lorentz_int
  real, dimension(:,:,:,:,:), allocatable ::  d_rbf_rs_int
  real, dimension(:,:,:), allocatable     ::  d_rbf_velint
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Entropy containers:
  !
  real, dimension(:,:), allocatable :: entropy
  !
  real, dimension(:,:,:,:), allocatable :: RHS_dr
  real, dimension(:,:,:,:), allocatable :: RHS_dt
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Species velocities
  !
  real, allocatable, dimension(:,:,:,:) :: v_para
  real, allocatable, dimension(:,:,:,:) :: v_perp
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Omegas:
  !
  real, allocatable, dimension(:,:,:,:) :: omega_d1
  real, allocatable, dimension(:,:,:,:) :: omega_star
  !
  real, allocatable, dimension(:,:,:,:) :: omega_dr
  !
  real, allocatable, dimension(:,:,:,:) :: v_theta
  real, allocatable, dimension(:) :: theta_int
  real, allocatable, dimension(:) :: theta_plot
  real, allocatable, dimension(:) :: theta_r0_plot
  !
  real, allocatable, dimension(:,:,:) :: lambda_b
  !---------------------------------------------------------

  !----------------------------------------------
  ! kx-ky spectrum
  !
  real, dimension(:), allocatable :: kxkyspec
  real, dimension(:), allocatable :: k_perp_squared
  !----------------------------------------------

  !------------------------------------------------
  ! Primitive fluxes:
  !
  real, dimension(:,:,:,:), allocatable :: nonlinear_flux_passing
  real, dimension(:,:,:,:), allocatable :: nonlinear_flux_trapped
  real, dimension(:,:), allocatable :: nonlinear_flux_momparts
  real, dimension(:,:,:,:,:), allocatable :: nonlinear_flux_velocity
  !
  ! Diffusivities:
  !
  real, dimension(:,:,:,:), allocatable :: diff_i
  real, dimension(:,:,:,:), allocatable :: diff_i_trapped
  real, dimension(:,:,:), allocatable :: diff
  real, dimension(:,:,:), allocatable :: diff_trapped
  real, dimension(:,:,:), allocatable :: diff_n
  real, dimension(:,:,:), allocatable :: s_diff_i
  real, dimension(:,:,:), allocatable :: s_diff_i_trapped
  real, dimension(:,:), allocatable :: s_diff
  real, dimension(:,:), allocatable :: s_diff_trapped
  real, dimension(:,:,:,:), allocatable :: sp_diff_i
  real, dimension(:,:,:,:), allocatable :: sp_diff_i_trapped
  real, dimension(:,:,:), allocatable :: sp_diff
  real, dimension(:,:,:), allocatable :: sp_diff_trapped
  !
  ! gyroBohm fluxes:
  !
  real, dimension(:,:,:,:), allocatable :: gbflux_i
  real, dimension(:,:,:,:), allocatable :: gbflux_i_trapped
  real, dimension(:,:,:), allocatable :: gbflux
  real, dimension(:,:), allocatable :: gbflux_mom
  real, dimension(:,:,:), allocatable :: gbflux_trapped
  real, dimension(:,:,:), allocatable :: gbflux_n
  !
  ! Additional diffusivities
  !
  real, dimension(1) :: phi_squared_QL_n
  real, dimension(3) :: g_squared_QL_n
  !
  real, dimension(:), allocatable :: Tr_p
  real, dimension(:), allocatable :: Eng_p
  !
  real, dimension(:,:,:,:), allocatable :: diff_vec
  real, dimension(:,:,:,:), allocatable :: gbflux_vec
  integer :: output_flag
  !------------------------------------------------

  !------------------------------------------------------
  ! Merged variables
  !
  real :: nu_source
  !
  real, dimension(:,:,:), allocatable :: h0_eq
  real, dimension(:,:,:), allocatable :: h0_mod
  !
  real, dimension(:,:), allocatable :: h0_n
  real, dimension(:,:), allocatable :: h0_e
  real, dimension(:,:), allocatable :: source_n
  real, dimension(:,:), allocatable :: source_e
  !
  real :: total_memory
  real, dimension(:), allocatable :: krho_collect
  !
  complex, dimension(:,:), allocatable :: z_gyro
  !
  integer, dimension(:), allocatable :: c_map
  !------------------------------------------------------

  !----------------------------------------------
  ! LAPACK
  !
  integer :: info
  integer, dimension(:), allocatable :: i_piv
  !-----------------------------------------------

  !------------------------------------------------------
  ! MPI
  !
  integer :: i_proc
  integer :: n_proc
  !
  integer :: i_proc_w
  !
  integer :: n_proc_1
  integer :: n_proc_2
  !
  integer :: i_proc_1
  integer :: i_proc_2
  !
  integer :: i_group_1
  integer :: i_group_2
  !
  integer :: GYRO_COMM_WORLD
  integer :: GKEIGEN_J_SUBSET
  integer :: NEW_COMM_1
  integer :: NEW_COMM_2
  integer :: MUMPS_COMM
  !
  integer :: i_err
  !
  integer :: msplit_SSUB
  integer :: nv1_SSUB
  integer :: nv2_SSUB
  !
  integer, dimension(:), allocatable :: recv_status
  !------------------------------------------------------

  !------------------------------------------------------
  ! UMFPACK control
  !
  real, dimension(4,10) :: cntl
  real, dimension(20) :: rinfo
  integer, dimension(4,20) :: icntl
  integer, dimension(40) :: uinfo
  integer, dimension(4,20) :: keep  
  !------------------------------------------------------

  !------------------------------------------
  ! Precision outut file
  !
  character (len=80) :: precfile='prec.out'
  !------------------------------------------

  !---------------------------------------------------------
  ! Path to INPUT, read in the get_inputpath subroutine
  !
  character(len=80) :: path
  !---------------------------------------------------------

end module gyro_globals
