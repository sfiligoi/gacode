module cgyro_globals

  use, intrinsic :: iso_c_binding

  !---------------------------------------------------------------
  ! Input parameters:
  !
  integer :: n_energy
  integer :: n_xi
  integer :: n_theta
  integer :: n_radial
  integer :: n_toroidal
  integer :: n_field
  integer :: e_max
  real    :: delta_t
  real    :: max_time
  integer :: print_step
  real    :: freq_tol
  integer :: restart_write
  integer :: restart_mode
  real    :: up_radial
  integer :: up_radial_n
  real    :: up_theta
  real    :: up_xi
  integer :: implicit_flag
  real    :: ky
  integer :: box_size
  integer :: silent_flag
  integer :: equilibrium_model
  integer :: collision_model
  integer :: collision_mom_restore
  integer :: collision_ene_restore
  integer :: collision_ene_diffusion
  integer :: collision_kperp
  integer :: collision_field_model
  integer :: collision_trap_model
  integer :: zf_test_flag 
  integer :: nonlinear_flag 
  integer :: nonlinear_method
  real :: te_ade
  real :: ne_ade
  real :: masse_ade
  real :: lambda_debye
  integer :: test_flag
  real :: amp
  real :: gamma_e
  !
  ! Geometry input
  !
  real :: rmin
  real :: rmaj
  real :: q
  real :: s
  real :: shift
  real :: kappa
  real :: s_kappa
  real :: delta
  real :: s_delta
  real :: zeta
  real :: s_zeta
  real :: zmag
  real :: s_zmag
  real :: beta_star
  real :: betae_unit
  !
  ! Species parameters
  !
  integer :: n_species
  real :: nu_ee_in
  integer, dimension(6) :: z
  real, dimension(6) :: mass
  real, dimension(6) :: dens
  real, dimension(6) :: temp
  real, dimension(6) :: dlnndr
  real, dimension(6) :: dlntdr
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! MPI variables and pointers
  ! 
  integer :: i_err
  integer :: i_proc
  integer :: i_proc_1
  integer :: i_proc_2
  integer :: n_proc
  integer :: n_proc_1
  integer :: n_proc_2
  integer :: i_group_1
  integer :: i_group_2
  integer :: CGYRO_COMM_WORLD
  integer :: NEW_COMM_1
  integer :: NEW_COMM_2
  integer :: nv1,nv2,nc1,nc2
  integer :: nsplit
  integer, dimension(:), allocatable :: recv_status
  !
  ! Pointers
  !
  integer :: nv,iv
  integer :: nv_loc,iv_loc
  integer :: nc,ic
  integer :: nc_loc,ic_loc
  integer, dimension(:), allocatable :: ie_v
  integer, dimension(:), allocatable :: ix_v
  integer, dimension(:), allocatable :: is_v
  integer, dimension(:), allocatable :: ir_c
  integer, dimension(:), allocatable :: it_c
  integer, dimension(:,:), allocatable :: ic_c
  integer, dimension(:,:,:), allocatable :: iv_v
  !
  integer :: n
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Constants
  !
  real, parameter    :: pi   = 3.1415926535897932
  complex, parameter :: i_c  = (0.0,1.0)
  real, parameter    :: num1 = 1.0
  real, parameter    :: num0 = 0.0
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! I/O and error management variables
  !
  character(len=80) :: path
  character(len=18) :: runfile_err     = 'out.cgyro.err'
  character(len=18) :: runfile_info    = 'out.cgyro.info'
  character(len=18) :: runfile_mpi     = 'out.cgyro.mpi'
  character(len=18) :: runfile_restart = 'out.cgyro.restart'
  character(len=18) :: runfile_hb      = 'out.cgyro.hb'
  character(len=18) :: runfile_caphb   = 'out.cgyro.caphb'
  character(len=18) :: runfile_grids   = 'out.cgyro.grids'
  character(len=18) :: runfile_time    = 'out.cgyro.time'
  character(len=18) :: runfile_freq    = 'out.cgyro.freq'
  character(len=14), dimension(3)  :: runfile_field = &
       (/'out.cgyro.phi ','out.cgyro.apar','out.cgyro.bpar'/)
  character(len=15), dimension(3)  :: runfile_fieldb = &
       (/'out.cgyro.phib ','out.cgyro.aparb','out.cgyro.bparb'/)
  character(len=16), dimension(2)  :: runfile_flux = &
       (/'out.cgyro.flux_n','out.cgyro.flux_e'/)
  character(len=18), dimension(3)  :: runfile_power = &
       (/'out.cgyro.pwr_phi ','out.cgyro.pwr_apar','out.cgyro.pwr_bpar'/)
  integer, parameter :: io=1
  !
  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message
  !
  integer :: io_control
  integer :: signal
  !
  ! Standard precision for IO 
  character(len=8) :: fmtstr='(es11.4)'
  character(len=14) :: fmtstr2='(2(es11.4,1x))'
  !----------------------------------------------------

  !---------------------------------------------------------------
  ! Time stepping
  integer :: i_time
  integer :: n_time
  complex :: freq
  complex :: freq_err
  real :: gtime
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Physics variables
  !
  integer :: ae_flag
  integer :: is_ele
  real :: dens_ele
  real :: temp_ele
  real :: mass_ele
  !
  real, dimension(6) :: vth  
  real, dimension(6) :: nu
  real :: rho
  real :: k_theta
  real :: length
  real :: omega_eb
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Numerical/work arrays and dimensions
  !
  ! Velocity space
  ! 
  integer, dimension(:), allocatable :: indx_xi, px
  real, dimension(:), allocatable :: energy, w_e
  real, dimension(:), allocatable :: xi, w_xi
  real, dimension(:,:), allocatable :: xi_deriv_mat, xi_lor_mat, xi_upderiv_mat
  real, dimension(:,:), allocatable :: e_deriv1_mat, e_deriv2_mat
  !
  ! Parallel streaming
  !
  real, dimension(:), allocatable :: theta
  real, dimension(-2:2) :: uderiv
  real, dimension(-2:2) :: cderiv
  integer, dimension(:), allocatable :: thcyc
  integer, dimension(:,:,:), allocatable :: rcyc
  integer, dimension(:), allocatable :: pcyc
  integer, dimension(:), allocatable :: ncyc
  complex, dimension(:,:,:), allocatable :: dtheta
  complex, dimension(:,:,:), allocatable :: dtheta_up
  !
  ! Distributions
  !
  complex, dimension(:,:,:), allocatable :: rhs
  complex, dimension(:,:), allocatable :: h_x
  complex, dimension(:,:), allocatable :: h0_x
  complex, dimension(:,:), allocatable :: psi
  complex, dimension(:,:,:), allocatable :: f_nl
  complex, dimension(:,:,:), allocatable :: g_nl
  complex, dimension(:,:), allocatable :: omega_cap_h
  complex, dimension(:,:), allocatable :: omega_h
  complex, dimension(:,:,:), allocatable :: omega_s
  complex, dimension(:,:), allocatable :: cap_h_c
  complex, dimension(:,:), allocatable :: cap_h_ct
  complex, dimension(:,:), allocatable :: cap_h_v
  complex, dimension(:,:), allocatable :: cap_h_v_prime
  real, dimension(:,:), allocatable :: j0_c
  real, dimension(:,:), allocatable :: j0_v
  !
  ! Fields
  !
  complex, dimension(:,:,:), allocatable :: field
  complex, dimension(:,:,:), allocatable :: field_loc
  complex, dimension(:,:,:), allocatable :: field_old
  complex, dimension(:,:,:), allocatable :: field_old2
  complex, dimension(:,:,:), allocatable :: field_old3
  !
  ! Nonlinear fluxes
  real, dimension(:,:), allocatable :: flux_loc
  real, dimension(:,:), allocatable :: flux
  real, dimension(:,:), allocatable :: power
  !
  type(C_PTR) :: plan_r2c
  type(C_PTR) :: plan_c2r
  !  
  ! Work arrays
  !
  complex, dimension(:,:), allocatable :: f_balloon
  real :: field_error
  !
  ! LAPACK work arrays 
  !
  real, dimension(:), allocatable :: work  
  integer, dimension(:), allocatable :: i_piv
  integer :: info
  !
  ! Implicit streaminggk/field matrices
  !
  complex, dimension(:,:,:), allocatable :: gkmat
  integer, dimension(:,:), allocatable   :: i_piv_gk
  complex, dimension(:,:), allocatable   :: gkvec
  complex, dimension(:,:), allocatable   :: fieldmat
  integer, dimension(:,:), allocatable   :: idfield
  integer, dimension(:),   allocatable   :: i_piv_field
  complex, dimension(:),   allocatable   :: fieldvec, fieldvec_loc
  ! umfpack
  integer, parameter :: gkmatsolve_flag=1
  real,    dimension(:,:), allocatable :: gksp_cntl
  integer, dimension(:,:), allocatable :: gksp_icntl, gksp_keep
  real,    dimension(20) ::  gksp_rinfo
  integer, dimension(40) ::  gksp_uinfo
  complex, dimension(:,:), allocatable :: gksp_mat
  integer, dimension(:,:), allocatable :: gksp_indx
  complex, dimension(:), allocatable   :: gksvec, gkwvec 
  integer :: gksp_nelem, gksp_nmax
  !
  !
  ! Some field solve parameters
  !
  real :: sum_den_h
  real, dimension(:,:), allocatable :: sum_den_x, sum_cur_x
  !
  ! n=0 test variables
  !
  real, dimension(:,:,:), allocatable :: hzf, xzf 
  real, dimension(:), allocatable :: pvec_in, pvec_outr, pvec_outi
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  integer :: geo_ny_in
  real, dimension(8,0:16) :: geo_yin_in
  !
  integer :: geo_numeq_flag
  integer :: geo_ny
  real, dimension(:,:), allocatable :: geo_yin
  !---------------------------------------------------------------

end module cgyro_globals
