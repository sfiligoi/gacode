module cgyro_globals

  !---------------------------------------------------------------
  ! local MPI variables
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
  !---------------------------------------------------------------

  real, parameter    :: pi   = 3.1415926535897932
  complex, parameter :: i_c  = (0.0,1.0)
  real, parameter    :: num1 = 1.0
  real, parameter    :: num0 = 0.0

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
  real    :: freq_tol
  integer :: restart_write
  integer :: restart_mode
  real    :: up_radial
  integer :: up_radial_n
  real    :: up_theta
  real    :: ky
  real    :: box_size
  integer :: silent_flag
  integer :: equilibrium_model
  integer :: collision_model 
  integer :: ae_flag
  integer :: zf_test_flag 
  real :: te_ade
  real :: ne_ade
  real :: lambda_debye

  ! Geometry input

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

  ! Species parameters

  integer :: n_species
  real :: nu_1_in
  integer, dimension(6) :: z
  real, dimension(6) :: mass
  real, dimension(6) :: dens
  real, dimension(6) :: temp
  real, dimension(6) :: dlnndr
  real, dimension(6) :: dlntdr

  !---------------------------------------------------------------

  real, dimension(6) :: vth  
  real, dimension(6) :: nu
  real :: rho
  real :: k_theta_rho
  real :: k_theta
  real :: r_length_rho
  real :: r_length_inv

  character(len=80)  :: runfile_restart = 'out.cgyro.restart'

  !---------------------------------------------------------------
  integer :: geo_ny_in
  real, dimension(8,0:16) :: geo_yin_in
  !
  integer :: geo_numeq_flag
  integer :: geo_ny
  real, dimension(:,:), allocatable :: geo_yin
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Time stepping
  integer :: itime, nt_step
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Path to input.cgyro, read in the get_inputpath subroutine
  character(len=80) :: path
  !---------------------------------------------------------------

  integer :: is_ele
  real :: dens_ele, temp_ele

  real, dimension(:), allocatable :: theta

  ! Numerical/work arrays and dimensions
  complex, dimension(:,:), allocatable :: h_x
  complex, dimension(:,:), allocatable :: cap_h_c
  complex, dimension(:,:), allocatable :: cap_h_ct
  complex, dimension(:,:), allocatable :: cap_h_v
  complex, dimension(:,:), allocatable :: cap_h_v_prime
  complex, dimension(:,:,:) , allocatable :: field
  complex, dimension(:,:,:) , allocatable :: field_loc
  complex, dimension(:,:,:) , allocatable :: field_old

  integer, dimension(:), allocatable :: indx_xi, indx_r
  real, dimension(:), allocatable :: energy, w_e
  real, dimension(:), allocatable :: xi, w_xi
  real, dimension(:,:), allocatable :: xi_deriv_mat, xi_lor_mat
  real, dimension(:,:), allocatable :: e_deriv1_mat, e_deriv2_mat

  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message

  ! output file
  integer :: output_flag=1
  character(len=80)  :: runfile = 'out.cgyro.run'
  integer :: io_cgyroout = 12

end module cgyro_globals
