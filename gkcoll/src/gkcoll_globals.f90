module gkcoll_globals

  !---------------------------------------------------------------
  ! local MPI variables
  ! 
  integer :: i_proc
  integer :: n_proc
  integer :: GKCOLL_COMM_WORLD
  !---------------------------------------------------------------

  real, parameter :: pi=3.1415926535897932

  !---------------------------------------------------------------
  ! Input parameters:
  !
  real :: rmin_in
  real :: rmaj_in
  real :: q_in
  real :: rho_in
  real :: shat_in
  real :: shift_in
  real :: kappa_in
  real :: s_kappa_in
  real :: delta_in
  real :: s_delta_in
  real :: zeta_in
  real :: s_zeta_in
  real :: zmag_in
  real :: s_zmag_in
  integer :: geo_ny_in
  real, dimension(8,0:16) :: geo_yin_in
  !
  real :: profile_delta_scale
  real :: profile_zeta_scale
  real :: profile_zmag_scale
  !
  integer :: ipccw_in
  integer :: btccw_in
  !
  real :: te_ade_in
  real :: ne_ade_in
  !
  integer, dimension(6) :: z_in 
  real, dimension(6) :: mass_in
  real, dimension(6) :: dens_in
  real, dimension(6) :: temp_in
  real, dimension(6) :: dlnndr_in
  real, dimension(6) :: dlntdr_in
  real, dimension(6) :: nu_in
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

  integer :: matsz_scalefac
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Models:
  !
  integer :: equilibrium_model
  integer :: collision_model 
  integer :: profile_model
  integer :: profile_temprescale_model
  integer :: profile_equilibrium_model
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
  ! Charge and mass
  !
  integer, dimension(:), allocatable :: z     ! charge (ns)
  real, dimension(:), allocatable :: mass     ! m/m_0  (ns)
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Path to INPUT, read in the get_inputpath subroutine
  character(len=80) :: path
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Profile functions
  !
  integer, parameter :: n_gr=1   ! radially-local only
  real, dimension(:), allocatable :: r
  real, dimension(:), allocatable :: rmaj
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: shat
  real, dimension(:), allocatable :: shift
  real, dimension(:), allocatable :: kappa
  real, dimension(:), allocatable :: s_kappa
  real, dimension(:), allocatable :: delta
  real, dimension(:), allocatable :: s_delta
  real, dimension(:), allocatable :: zeta
  real, dimension(:), allocatable :: s_zeta
  real, dimension(:), allocatable :: zmag
  real, dimension(:), allocatable :: s_zmag
  integer :: geo_numeq_flag
  integer :: geo_ny
  real, dimension(:,:,:), allocatable :: geo_yin
  !
  real, dimension(:), allocatable :: te_ade, ne_ade
  ! ele temp and dens -- used only if ade (for Poisson calculation)
  !
  ! (species-dependent)
  !
  real, dimension(:,:), allocatable :: dens   ! n/n_0  (ns,nr)
  real, dimension(:,:), allocatable :: temp   ! T/T_0  (ns,nr)
  real, dimension(:,:), allocatable :: vth    ! vth/vth_0=sqrt(temp/mass) 
  real, dimension(:,:), allocatable :: dlnndr ! a/Ln   (ns,nr)
  real, dimension(:,:), allocatable :: dlntdr ! a/LT   (ns,nr)
  real, dimension(:,:), allocatable :: nu     ! nu / (vth_0/a) (ns,nr)
  !

  !---------------------------------------------------------------

  real, dimension(:), allocatable :: theta

  ! Numerical/work arrays and dimensions
  real, dimension(:,:,:,:,:), allocatable :: f
  real, dimension(:,:) , allocatable :: phi

  real, dimension(:), allocatable :: indx_xi, indx_r

  ! normalizations for experimental profiles
  real  :: temp_norm_fac, charge_norm_fac
  real  :: a_meters
  real, dimension(:), allocatable :: dens_norm
  real, dimension(:), allocatable :: temp_norm
  real, dimension(:), allocatable :: vth_norm
  real, dimension(:), allocatable :: b_unit
  real, parameter :: mass_deuterium = 3.3452   ! (x 10-27 kg)

  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message
	
  ! output file
  character(len=80)  :: runfile_gkcollout = 'out.gkcoll.run'
  integer :: io_gkcollout = 12


end module gkcoll_globals
