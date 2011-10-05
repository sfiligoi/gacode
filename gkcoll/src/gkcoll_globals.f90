module gkcoll_globals

  !---------------------------------------------------------------
  ! local MPI variables
  ! 
  integer :: i_proc
  integer :: n_proc
  integer :: GKCOLL_COMM_WORLD
  !---------------------------------------------------------------

  real, parameter :: pi=3.1415926535897932
  complex, parameter :: i_c = (0.0,1.0)

  !---------------------------------------------------------------
  ! Input parameters:
  !
  real :: rmin
  real :: rmaj
  real :: q
  real :: shat
  real :: shift
  real :: kappa
  real :: s_kappa
  real :: delta
  real :: s_delta
  real :: zeta
  real :: s_zeta
  real :: zmag
  real :: s_zmag
  
  integer :: toroidal_model
  integer :: toroidal_num
  real :: rho
  real :: k_theta_rho
  real :: k_theta
  real :: r_length_inv

  !---------------------------------------------------------------
  integer :: geo_ny_in
  real, dimension(8,0:16) :: geo_yin_in
  !
  integer :: geo_numeq_flag
  integer :: geo_ny
  real, dimension(:,:), allocatable :: geo_yin
  !---------------------------------------------------------------

  !
  integer :: ipccw
  integer :: btccw
  !
  real :: te_ade
  real :: ne_ade
  !
  real :: lambda_debye
  real :: profile_lambda_debye_scale
  !
  integer, dimension(6) :: z
  real, dimension(6) :: mass
  real, dimension(6) :: dens
  real, dimension(6) :: temp
  real, dimension(6) :: dlnndr
  real, dimension(6) :: dlntdr
  real, dimension(6) :: nu
  ! 
  real, dimension(6) :: vth  
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Grid dimensions:
  !
  integer :: n_species
  integer :: n_energy
  integer :: n_xi
  integer :: n_theta
  integer :: n_radial
  real    :: e_max
  real    :: delta_t
  real    :: max_time
  real    :: freq_tol
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Models:
  !
  integer :: equilibrium_model
  integer :: collision_model 
  integer :: profile_model
  integer :: adiabatic_ele_model
  real    :: sign_q
  real    :: sign_bunit
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Output mode:
  !
  integer :: silent_flag
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Path to INPUT, read in the get_inputpath subroutine
  character(len=80) :: path
  !---------------------------------------------------------------

  integer :: is_ele
  real :: dens_ele, temp_ele

  real, dimension(:), allocatable :: theta

  ! Numerical/work arrays and dimensions
  complex, dimension(:,:,:,:,:), allocatable :: h_x
  complex, dimension(:,:,:,:,:), allocatable :: cap_h_x
  complex, dimension(:,:,:,:,:), allocatable :: cap_h_p
  complex, dimension(:,:) , allocatable :: phi
  complex, dimension(:,:) , allocatable :: phi_old

  integer, dimension(:), allocatable :: indx_xi, indx_r
  real, dimension(:), allocatable :: energy, w_e
  real, dimension(:), allocatable :: xi, w_xi

  ! normalizations for experimental profiles
  real  :: temp_norm_fac, charge_norm_fac
  real  :: a_norm
  real  :: dens_norm
  real  :: temp_norm
  real  :: vth_norm
  real  :: b_norm
  real, parameter :: mass_deuterium = 3.3452   ! (x 10-27 kg)

  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message
	
  ! output file
  character(len=80)  :: runfile_gkcollout = 'out.gkcoll.run'
  integer :: io_gkcollout = 12


end module gkcoll_globals
