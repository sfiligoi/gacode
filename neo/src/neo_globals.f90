module neo_globals

  !---------------------------------------------------------------
  ! local MPI variables
  ! 
  integer :: i_proc
  integer :: n_proc
  integer :: NEO_COMM_WORLD
  !---------------------------------------------------------------

  real, parameter :: pi=3.1415926535897932

  !---------------------------------------------------------------
  ! Input parameters:
  !
  real :: rmin_1_in
  real :: rmin_2_in
  real :: rmaj_in
  real :: dphi0dr_in
  real :: epar0_in
  real :: omega_rot_in
  real :: omega_rot_deriv_in
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

  real :: profile_dlnndr_1_scale
  real :: profile_dlnndr_2_scale
  real :: profile_dlnndr_3_scale
  real :: profile_dlnndr_4_scale
  real :: profile_dlnndr_5_scale
  real :: profile_dlnndr_6_scale

  real :: profile_dlntdr_1_scale
  real :: profile_dlntdr_2_scale
  real :: profile_dlntdr_3_scale
  real :: profile_dlntdr_4_scale
  real :: profile_dlntdr_5_scale
  real :: profile_dlntdr_6_scale
  real, dimension(6) :: dlnndr_scale
  real, dimension(6) :: dlntdr_scale
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
  !
  real :: nu_1_in
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Grid dimensions:
  !
  integer :: n_species
  integer :: n_energy
  integer :: n_xi
  integer :: n_theta
  integer :: n_radial

  integer :: matsz_scalefac
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Models:
  !
  integer :: sim_model
  integer :: equilibrium_model
  integer :: collision_model 
  integer :: profile_model
  integer :: profile_erad0_model
  integer :: profile_equilibrium_model
  integer :: rotation_model
  integer :: adiabatic_ele_model
  integer :: spitzer_model
  real    :: epar0_spitzer
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
  real, dimension(:), allocatable :: dphi0dr  ! dphi0/dr (e/T_0)*a (nr) 
  ! erad0 = -dphi0/dr * |grad r(theta=0)|
  real, dimension(:), allocatable :: epar0    ! inductive Efield (e/T_0)*a (nr)
  real, dimension(:), allocatable :: omega_rot  ! omega_rot / (vth_0/a) (nr)
  real, dimension(:), allocatable :: omega_rot_deriv ! domega_rot/d(r/a) (nr)

  !---------------------------------------------------------------

  real, dimension(:), allocatable :: theta

  ! Numerical/work arrays and dimensions

  integer :: n_row
  integer :: n_max

  real, dimension(:), allocatable :: g
  real, dimension(:), allocatable :: amat
  integer, dimension(:), allocatable :: amat_indx

  integer, dimension(:), allocatable :: is_indx, ie_indx, ix_indx, it_indx

  real, dimension(:,:), allocatable :: driftx, &
       driftxrot1, driftxrot2, driftxrot3
  ! driftx   = vdrift dot grad r

  ! normalizations for experimental profiles
  real  :: temp_norm_fac, charge_norm_fac
  real  :: a_meters
  real, dimension(:), allocatable :: dens_norm
  real, dimension(:), allocatable :: temp_norm
  real, dimension(:), allocatable :: vth_norm
  real, dimension(:), allocatable :: b_unit
  real, parameter :: mass_deuterium = 3.3452   ! (x 10-27 kg)

  ! theta derivative variables
  integer, dimension(-2:2) :: cderiv
  integer, dimension(:), allocatable :: thcyc

  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message
	
  ! output file
  character(len=80)  :: runfile_neoout = 'out.neo.run'
  integer :: io_neoout = 12

  ! output vectors

  ! (n_species_max, transport coeff)
  ! transport coeff: 1-> gamma, 2-> Q, 3->Pi, 4-> Q-omega*Pi
  !                  5-> vpol,  6-> vtor
  real, dimension(6,6) :: neo_dke_out=0.0

  ! species-independent transport coeff (currently just jpar)
  real                 :: neo_dke_1d_out= 0.0

  ! (n_species_max, transport coeff)  
  ! gyro-viscous fluxes: 1-> gamma, 2-> Q, 3->Pi, 4-> Q-omega*Pi
  real, dimension(6,4) :: neo_gv_out=0.0

  ! pure plasma theory
  ! Gamma_HH, Qi_HH, Qe_HH, Qi_CH, jpar_HH, jpar_S
  real, dimension(6)   :: neo_th_out=0.0

  ! (n_species_max, transport coeff)
  ! multi-species theory: 1-> gamma_HS, 2-> Q_HS
  real, dimension(6,2)   :: neo_thHS_out=0.0

end module neo_globals
