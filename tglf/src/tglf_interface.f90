!-------------------------------------------------------------------------
! tglf_interface.f90
!
! PURPOSE:
!  Provides interface description for TGLF.
!
! CALLING SEQUENCE:
!  call tglf_init(...)
!  set tglf_*_in variables
!  call tglf_run(...)
!  get tglf_*_out variables
!
! NOTES:
!  tglf_dump_flag:   TRUE = call dump routines, FALSE = do NOT call dump routines 
!  tglf_dump_suffix: appended to name of dump file 
!  tglf_dump_path:   specifies directory where dump files should be 
!
!  tglf_dump_local and tglf_dump_global are called by tglf_run
!
! NAMING CONVENTION:
!  Interface parameter names are formed from TGLF variable names by 
!  prepending 'tglf_' and appending '_in'
!
! INITIALIZATION:
!  We initialize input arrays/vectors to zero unlike original TGLF module
!
!-------------------------------------------------------------------------

module tglf_interface

  implicit none

  ! PRIVATE PARAMETERS
  ! Number of species (should correspond to nsm in tglf_dimensions)
  integer, parameter :: nsm = 6

  ! CONTROL PARAMETERS
  integer              :: tglf_geometry_flag_in = 1
  logical              :: tglf_dump_flag_in     = .false.
  character (len=64)   :: tglf_dump_suffix_in   = ''
  character (len=4096) :: tglf_dump_path_in     = ''

  ! INPUT PARAMETERS

  ! Data passed to: put_rare_switches
  real    :: tglf_alpha_quench_in   = 0.0
  real    :: tglf_park_in           = 1.0
  real    :: tglf_ghat_in           = 1.0
  real    :: tglf_gchat_in          = 1.0
  real    :: tglf_wd_zero_in        = 0.1
  real    :: tglf_linsker_factor_in = 0.0
  real    :: tglf_gradB_factor_in   = 0.0
  real    :: tglf_x_psi_in          = 0.0
  integer :: tglf_sat_rule_in       = 1 
  real    :: tglf_alpha_kx0_in      = 1.0

  ! Data passed to: put_switches
  logical :: tglf_iflux_in          = .true.
  logical :: tglf_use_bper_in       = .false.
  logical :: tglf_use_bpar_in       = .false.
  logical :: tglf_use_bisection_in  = .false.
  integer :: tglf_ibranch_in        = -1
  integer :: tglf_nmodes_in         = 2
  integer :: tglf_nbasis_max_in     = 4
  integer :: tglf_nbasis_min_in     = 1
  integer :: tglf_nxgrid_in         = 16
  integer :: tglf_nky_in            = 12
  integer :: tglf_kygrid_model_in   = 1
  integer :: tglf_xnu_model_in      = 2

  ! Data passed to: put_model_parameters
  logical :: tglf_adiabatic_elec_in = .false.
  real    :: tglf_alpha_p_in        = 0.0
  real    :: tglf_alpha_e_in        = 0.18
  real    :: tglf_theta_trapped_in  = 0.7
  real    :: tglf_xnu_factor_in     = 1.0
  real    :: tglf_debye_factor_in   = 1.0
  real    :: tglf_etg_factor_in     = 1.25
  real    :: tglf_filter_in         = 2.0

  ! Data passed to: put_species
  integer :: tglf_ns_in             = 2
  real    :: tglf_mass_in(nsm)      = 0.0
  real    :: tglf_zs_in(nsm)        = 0.0

  ! Data passed to: put_kys
  real    :: tglf_ky_in             = 0.3

  ! Data passed to: put_gaussian_width
  real    :: tglf_width_in          = 1.65
  real    :: tglf_width_min_in      = 0.3
  integer :: tglf_nwidth_in         = 21
  logical :: tglf_find_width_in     = .false.

  ! Data passed to: put_gradients
  real    :: tglf_rlns_in(nsm)      = 0.0
  real    :: tglf_rlts_in(nsm)      = 0.0
  real    :: tglf_vexb_shear_in     = 0.0
  real    :: tglf_vpar_shear_in     = 0.0

  ! Data passed to: put_averages
  real    :: tglf_taus_in(nsm)      = 0.0
  real    :: tglf_as_in(nsm)        = 0.0
  real    :: tglf_betae_in          = 0.0
  real    :: tglf_xnuei_in          = 0.0
  real    :: tglf_zeff_in           = 1.0
  real    :: tglf_debye_in          = 0.0
  real    :: tglf_vexb_mach_in      = 0.0
  real    :: tglf_vpar_mach_in      = 0.0

  ! Data passed to: put_eikonal
  logical :: tglf_new_eikonal_in    = .true.

  ! Data passed to: put_s_alpha_geometry
  real    :: tglf_rmin_sa_in        = 0.5
  real    :: tglf_rmaj_sa_in        = 3.0
  real    :: tglf_q_sa_in           = 2.0
  real    :: tglf_shat_sa_in        = 1.0
  real    :: tglf_alpha_sa_in       = 0.0
  real    :: tglf_xwell_sa_in       = 0.0
  real    :: tglf_theta0_sa_in      = 0.0
  integer :: tglf_b_model_sa_in     = 1
  integer :: tglf_ft_model_sa_in    = 1

  ! Data passed to: put_Miller_geometry
  real    :: tglf_rmin_loc_in       = 0.5
  real    :: tglf_rmaj_loc_in       = 3.0
  real    :: tglf_q_loc_in          = 2.0
  real    :: tglf_q_prime_loc_in    = 16.0
  real    :: tglf_p_prime_loc_in    = 0.0
  real    :: tglf_shift_loc_in      = 0.0
  real    :: tglf_kappa_loc_in      = 1.0
  real    :: tglf_s_kappa_loc_in    = 0.0
  real    :: tglf_delta_loc_in      = 0.0
  real    :: tglf_s_delta_loc_in    = 0.0

  ! OUTPUT PARAMETERS
  real :: tglf_elec_pflux_out = 0.0
  real :: tglf_elec_eflux_out = 0.0
  real :: tglf_elec_eflux_low_out = 0.0
  real :: tglf_elec_mflux_out = 0.0

  real, dimension(5) :: tglf_ion_pflux_out = 0.0
  real, dimension(5) :: tglf_ion_eflux_out = 0.0
  real, dimension(5) :: tglf_ion_eflux_low_out = 0.0
  real, dimension(5) :: tglf_ion_mflux_out = 0.0

contains

  ! Dump LOCAL INTERFACE variables
  subroutine tglf_dump_local()

    implicit none

    integer :: ioerr

    open(unit=1,file=trim(tglf_dump_path_in)//'tglf_local.dump'//trim(tglf_dump_suffix_in),&
         status='replace',iostat=ioerr)

    write(1,*) 'new_eikonal = ', tglf_new_eikonal_in
    write(1,*) 'find_width = ', tglf_find_width_in
    write(1,*) 'nwidth = ', tglf_nwidth_in
    write(1,*) 'width = ', tglf_width_in
    write(1,*) 'width_min = ', tglf_width_min_in
    write(1,*) 'ky = ', tglf_ky_in
    write(1,*) 'nky = ', tglf_nky_in
    write(1,*) 'nmodes = ', tglf_nmodes_in
    write(1,*) 'ns = ', tglf_ns_in
    write(1,*) 'mass = ', tglf_mass_in(:)
    write(1,*) 'zs = ', tglf_zs_in(:)
    write(1,*) 'iflux = ', tglf_iflux_in
    write(1,*) 'use_bper = ', tglf_use_bper_in
    write(1,*) 'use_bpar = ', tglf_use_bpar_in
    write(1,*) 'use_bisection = ', tglf_use_bisection_in
    write(1,*) 'ibranch = ', tglf_ibranch_in
    write(1,*) 'nbasis_max = ', tglf_nbasis_max_in
    write(1,*) 'nbasis_min = ', tglf_nbasis_min_in
    write(1,*) 'nxgrid = ', tglf_nxgrid_in
    write(1,*) 'park = ', tglf_park_in
    write(1,*) 'ghat = ', tglf_ghat_in
    write(1,*) 'gchat = ', tglf_gchat_in
    write(1,*) 'wd_zero = ', tglf_wd_zero_in
    write(1,*) 'linsker_factor = ', tglf_linsker_factor_in
    write(1,*) 'gradB_factor = ', tglf_gradB_factor_in
    write(1,*) 'filter = ', tglf_filter_in
    write(1,*) 'x_psi = ', tglf_x_psi_in
    write(1,*) 'adiabatic_elec = ', tglf_adiabatic_elec_in
    write(1,*) 'alpha_p = ', tglf_alpha_p_in
    write(1,*) 'alpha_e = ', tglf_alpha_e_in
    write(1,*) 'theta_trapped = ', tglf_theta_trapped_in
    write(1,*) 'xnu_factor = ', tglf_xnu_factor_in
    write(1,*) 'debye_factor = ', tglf_debye_factor_in
    write(1,*) 'etg_factor = ', tglf_etg_factor_in
    write(1,*) 'rlns = ', tglf_rlns_in(:)
    write(1,*) 'rlts = ', tglf_rlts_in(:)
    write(1,*) 'vexb_shear = ', tglf_vexb_shear_in
    write(1,*) 'vpar_shear = ', tglf_vpar_shear_in
    write(1,*) 'as = ', tglf_as_in(:)
    write(1,*) 'taus = ', tglf_taus_in(:)
    write(1,*) 'betae = ', tglf_betae_in
    write(1,*) 'xnuei = ', tglf_xnuei_in
    write(1,*) 'zeff = ', tglf_zeff_in
    write(1,*) 'debye = ', tglf_debye_in

    if (tglf_geometry_flag_in == 1) then

       write(1,*) 'rmin_loc = ', tglf_rmin_loc_in
       write(1,*) 'rmaj_loc = ', tglf_rmaj_loc_in
       write(1,*) 'q_loc = ', tglf_q_loc_in
       write(1,*) 'shift_loc = ', tglf_shift_loc_in
       write(1,*) 'kappa_loc = ', tglf_kappa_loc_in
       write(1,*) 's_kappa_loc = ', tglf_s_kappa_loc_in
       write(1,*) 'delta_loc = ', tglf_delta_loc_in
       write(1,*) 's_delta_loc = ', tglf_s_delta_loc_in
       write(1,*) 'p_prime_loc = ', tglf_p_prime_loc_in
       write(1,*) 'q_prime_loc = ', tglf_q_prime_loc_in

    else

       write(1,*) 'rmin_sa = ', tglf_rmin_sa_in
       write(1,*) 'rmaj_sa = ', tglf_rmaj_sa_in
       write(1,*) 'q_sa = ', tglf_q_sa_in
       write(1,*) 'shat_sa = ', tglf_shat_sa_in
       write(1,*) 'alpha_sa = ', tglf_alpha_sa_in
       write(1,*) 'xwell_sa = ', tglf_xwell_sa_in
       write(1,*) 'theta0_sa = ', tglf_theta0_sa_in
       write(1,*) 'b_model_sa = ', tglf_b_model_sa_in
       write(1,*) 'ft_model_sa = ', tglf_ft_model_sa_in

    endif

    close(1)

  end subroutine tglf_dump_local

  ! Dump GLOBAL MODULE variables
  subroutine tglf_dump_global()

    use tglf_internal_interface

    implicit none

    integer :: ioerr

    open(unit=1,file=trim(tglf_dump_path_in)//'tglf_global.dump'//trim(tglf_dump_suffix_in),&
         status='replace',iostat=ioerr)

    write(1,*) 'new_eikonal = ', new_eikonal_in
    write(1,*) 'find_width = ', find_width_in
    write(1,*) 'nwidth = ', nwidth_in
    write(1,*) 'width = ', width_in
    write(1,*) 'width_min = ', width_min_in
    write(1,*) 'ky = ', ky_in
    write(1,*) 'nky = ', nky_in
    write(1,*) 'nmodes = ', nmodes_in
    write(1,*) 'ns = ', ns_in
    write(1,*) 'mass = ', mass_in(:)
    write(1,*) 'zs = ', zs_in(:)
    write(1,*) 'iflux = ', iflux_in
    write(1,*) 'use_bper = ', use_bper_in
    write(1,*) 'use_bpar = ', use_bpar_in
    write(1,*) 'use_bisection = ', use_bisection_in
    write(1,*) 'ibranch = ', ibranch_in
    write(1,*) 'nbasis_max = ', nbasis_max_in
    write(1,*) 'nbasis_min = ', nbasis_min_in
    write(1,*) 'nxgrid = ', nxgrid_in
    write(1,*) 'park = ', park_in
    write(1,*) 'ghat = ', ghat_in
    write(1,*) 'gchat = ', gchat_in
    write(1,*) 'wd_zero = ', wd_zero_in
    write(1,*) 'linsker_factor = ', linsker_factor_in
    write(1,*) 'gradB_factor = ', gradB_factor_in
    write(1,*) 'filter = ', filter_in
    write(1,*) 'x_psi = ', x_psi_in
    write(1,*) 'adiabatic_elec = ', adiabatic_elec_in
    write(1,*) 'alpha_p = ', alpha_p_in
    write(1,*) 'alpha_e = ', alpha_e_in
    write(1,*) 'theta_trapped = ', theta_trapped_in
    write(1,*) 'xnu_factor = ', xnu_factor_in
    write(1,*) 'debye_factor = ', debye_factor_in
    write(1,*) 'etg_factor = ', etg_factor_in
    write(1,*) 'rlns = ', rlns_in(:)
    write(1,*) 'rlts = ', rlts_in(:)
    write(1,*) 'vexb_shear = ', vexb_shear_in
    write(1,*) 'vpar_shear = ', vpar_shear_in
    write(1,*) 'as = ', as_in(:)
    write(1,*) 'taus = ', taus_in(:)
    write(1,*) 'betae = ', betae_in
    write(1,*) 'xnuei = ', xnuei_in
    write(1,*) 'zeff = ', zeff_in
    write(1,*) 'debye = ', debye_in

    if (tglf_geometry_flag_in == 1 ) then

       write(1,*) 'rmin_loc = ', rmin_loc
       write(1,*) 'rmaj_loc = ', rmaj_loc
       write(1,*) 'q_loc = ', q_loc
       write(1,*) 'shift_loc = ', shift_loc
       write(1,*) 'kappa_loc = ', kappa_loc
       write(1,*) 's_kappa_loc = ', s_kappa_loc
       write(1,*) 'delta_loc = ', delta_loc
       write(1,*) 's_delta_loc = ', s_delta_loc
       write(1,*) 'p_prime_loc = ', p_prime_loc
       write(1,*) 'q_prime_loc = ', q_prime_loc

    else

       write(1,*) 'rmin_sa = ', rmin_sa
       write(1,*) 'rmaj_sa = ', rmaj_sa
       write(1,*) 'q_sa = ', q_sa
       write(1,*) 'shat_sa = ', shat_sa
       write(1,*) 'alpha_sa = ', alpha_sa
       write(1,*) 'xwell_sa = ', xwell_sa
       write(1,*) 'theta0_sa = ', theta0_sa
       write(1,*) 'b_model_sa = ', b_model_sa
       write(1,*) 'ft_model_sa = ', ft_model_sa

    endif

    close(1)

  end subroutine tglf_dump_global

end module tglf_interface
