!-------------------------------------------------------------------------
! gftm_interface.f90
!
! PURPOSE:
!  Provides interface description for gftm.
!
! CALLING SEQUENCE:
!  call gftm_init(...)
!  set gftm_*_in variables
!  call gftm_run(...)
!  get gftm_*_out variables
!
! NOTES:
!  gftm_path:        working directory
!  gftm_dump_flag:   TRUE = call dump routines, FALSE = do NOT call dump routines
!  gftm_dump_suffix: appended to name of dump file
!
!  gftm_dump_local and gftm_dump_global are called by gftm_run
!
! NAMING CONVENTION:
!  Interface parameter names are formed from gftm variable names by
!  prepending 'gftm_' and appending '_in'
!
! INITIALIZATION:
!  We initialize input arrays/vectors to zero unlike original gftm module
!
!-------------------------------------------------------------------------

module gftm_interface

  use gftm_max_dimensions

  !  NOTE: gftm_max_dimensions defines
  !   nsm=6
  !   maxmodes=4

  implicit none

  ! CONTROL PARAMETERS
  character (len=256)  :: gftm_path_in       = ''
  character (len=20)   :: file_dump_local = 'out.gftm.localdump'
  logical              :: gftm_dump_flag_in  = .false.
  logical              :: gftm_quiet_flag_in = .true.
  integer              :: gftm_test_flag_in  = 0

! units switch
  character (len=8) :: gftm_units_in = 'GYRO'

! INPUT PARAMETERS
  logical :: gftm_use_transport_model_in = .true.
  integer :: gftm_geometry_flag_in = 1
  integer :: gftm_write_wavefunction_flag_in=0

  ! Data passed to: put_signs
  real    :: gftm_sign_bt_in        = 1.0
  real    :: gftm_sign_it_in        = 1.0

  ! Data passed to: put_rare_switches
  real    :: gftm_theta_trapped_in  = 0.7
  real    :: gftm_wdia_trapped_in   = 0.0
  real    :: gftm_park_in           = 1.0
  real    :: gftm_ghat_in           = 1.0
  real    :: gftm_gchat_in          = 1.0
  real    :: gftm_wd_zero_in        = 0.1
  real    :: gftm_linsker_factor_in = 0.0
  real    :: gftm_gradb_factor_in   = 0.0
  real    :: gftm_filter_in         = 2.0
  real    :: gftm_damp_psi_in       = 0.0
  real    :: gftm_damp_sig_in       = 0.0

  ! Data passed to: put_switches
  logical :: gftm_iflux_in          = .true.
  logical :: gftm_use_bper_in       = .false.
  logical :: gftm_use_bpar_in       = .false.
  logical :: gftm_use_mhd_rule_in   = .true.
  logical :: gftm_use_bisection_in  = .true.
  logical :: gftm_use_inboard_detrapped_in = .false.
  logical :: gftm_use_ave_ion_grid_in = .false.
  integer :: gftm_ibranch_in        = -1
  integer :: gftm_nmodes_in         = 2
  integer :: gftm_nbasis_max_in     = 4
  integer :: gftm_nbasis_min_in     = 2
  integer :: gftm_nxgrid_in         = 16
  integer :: gftm_nky_in            = 12

  ! Data passed to: put_model_parameters
  logical :: gftm_adiabatic_elec_in = .false.
  real    :: gftm_alpha_mach_in     = 0.0
  real    :: gftm_alpha_e_in        = 1.0
  real    :: gftm_alpha_p_in        = 1.0
  real    :: gftm_alpha_quench_in   = 0.0
  real    :: gftm_alpha_zf_in       = 1.0
  real    :: gftm_xnu_factor_in     = 1.0
  real    :: gftm_debye_factor_in   = 1.0
  real    :: gftm_etg_factor_in     = 1.25
  real    :: gftm_rlnp_cutoff_in     = 18.0
  integer :: gftm_sat_rule_in       = 0
  integer :: gftm_kygrid_model_in   = 1
  integer :: gftm_xnu_model_in      = 2
  integer :: gftm_vpar_model_in     = 0
  integer :: gftm_vpar_shear_model_in = 1

  ! Data passed to: put_species
  integer :: gftm_ns_in             = 2
  real    :: gftm_mass_in(nsm)      = 0.0
  real    :: gftm_zs_in(nsm)        = 0.0

  ! Data passed to: put_kys
  real    :: gftm_ky_in             = 0.3

  ! Data passed to: put_gaussian_width
  real    :: gftm_width_in          = 1.65
  real    :: gftm_width_min_in      = 0.3
  integer :: gftm_nwidth_in         = 21
  logical :: gftm_find_width_in     = .true.

  ! Data passed to: put_gradients
  real    :: gftm_rlns_in(nsm)      = 0.0
  real    :: gftm_rlts_in(nsm)      = 0.0
  real    :: gftm_vpar_shear_in(nsm)= 0.0
  real    :: gftm_vexb_shear_in     = 0.0

  ! Data passed to: put_profile_shear
  real    :: gftm_vns_shear_in(nsm) = 0.0
  real    :: gftm_vts_shear_in(nsm) = 0.0

  ! Data passed to: put_averages
  real    :: gftm_taus_in(nsm)    = 0.0
  real    :: gftm_as_in(nsm)      = 0.0
  real    :: gftm_vpar_in(nsm)    = 0.0
  real    :: gftm_vexb_in         = 0.0
  real    :: gftm_betae_in        = 0.0
  real    :: gftm_xnue_in         = 0.0
  real    :: gftm_zeff_in         = 1.0
  real    :: gftm_debye_in        = 0.0

  ! Data passed to: put_eikonal
  logical :: gftm_new_eikonal_in    = .true.

  ! Data passed to: put_s_alpha_geometry
  real    :: gftm_rmin_sa_in        = 0.5
  real    :: gftm_rmaj_sa_in        = 3.0
  real    :: gftm_q_sa_in           = 2.0
  real    :: gftm_shat_sa_in        = 1.0
  real    :: gftm_alpha_sa_in       = 0.0
  real    :: gftm_xwell_sa_in       = 0.0
  real    :: gftm_theta0_sa_in      = 0.0
  integer :: gftm_b_model_sa_in     = 1
  integer :: gftm_ft_model_sa_in    = 1

  ! Data passed to: put_Miller_geometry
  real    :: gftm_rmin_loc_in       = 0.5
  real    :: gftm_rmaj_loc_in       = 3.0
  real    :: gftm_zmaj_loc_in       = 0.0
  real    :: gftm_drmindx_loc_in    = 1.0
  real    :: gftm_drmajdx_loc_in    = 0.0
  real    :: gftm_dzmajdx_loc_in    = 0.0
  real    :: gftm_kappa_loc_in      = 1.0
  real    :: gftm_s_kappa_loc_in    = 0.0
  real    :: gftm_delta_loc_in      = 0.0
  real    :: gftm_s_delta_loc_in    = 0.0
  real    :: gftm_zeta_loc_in       = 0.0
  real    :: gftm_s_zeta_loc_in     = 0.0
  real    :: gftm_q_loc_in          = 2.0
  real    :: gftm_q_prime_loc_in    = 16.0
  real    :: gftm_p_prime_loc_in    = 0.0
  real    :: gftm_kx0_loc_in        = 0.0

  ! Data passed to Fourier_geometry
  integer :: gftm_nfourier_in        = 16
  real :: gftm_q_fourier_in          = 2.0
  real :: gftm_q_prime_fourier_in    = 16.0
  real :: gftm_p_prime_fourier_in    = 0.0
  real :: gftm_fourier_in(8,0:max_fourier)=0.0

  ! Data passed to: put_ELITE_geometry
  integer  :: gftm_n_elite_in   = 700
  real     :: gftm_q_elite_in   = 2.0
  real     :: gftm_q_prime_elite_in = 16.0
  real     :: gftm_R_elite_in(max_ELITE)
  real     :: gftm_Z_elite_in(max_ELITE)
  real     :: gftm_Bp_elite_in(max_ELITE)
  
  ! TRANSPORT OUTPUT PARAMETERS
  real :: gftm_elec_pflux_out = 0.0
  real :: gftm_elec_eflux_out = 0.0
  real :: gftm_elec_eflux_low_out = 0.0
  real :: gftm_elec_mflux_out = 0.0
  real :: gftm_elec_expwd_out = 0.0

  real, dimension(nsm-1) :: gftm_ion_pflux_out = 0.0
  real, dimension(nsm-1) :: gftm_ion_eflux_out = 0.0
  real, dimension(nsm-1) :: gftm_ion_eflux_low_out = 0.0
  real, dimension(nsm-1) :: gftm_ion_mflux_out = 0.0
  real, dimension(nsm-1) :: gftm_ion_expwd_out = 0.0

  real, dimension(nsm-1, 3) :: gftm_particle_flux_out = 0.0
  real, dimension(nsm-1, 3) :: gftm_energy_flux_out = 0.0
  real, dimension(nsm-1, 3) :: gftm_stress_tor_out = 0.0
  real, dimension(nsm-1, 3) :: gftm_stress_par_out = 0.0
  real, dimension(nsm-1, 3) :: gftm_exchange_out = 0.0
  
  ! LINEAR OUTPUT PARAMETERS
  complex :: gftm_eigenvalue_out(maxmodes)

  ! GYRO gftm input
  real, allocatable, dimension(:, :, :) :: gftm_eigenvalue_spectrum_out
  real, allocatable, dimension(:) :: gftm_ky_spectrum_out, gftm_dky_spectrum_out

  real, allocatable, dimension(:, :, :) :: gftm_field_spectrum_out
  real, allocatable, dimension(:, :, :, :, :) :: gftm_flux_spectrum_out
  
  
  
  ! DIAGNOSTIC OUTPUT PARAMETERS
  real :: interchange_DR = 0.0
  real :: interchange_DM = 0.0

  ! ERROR OUTPUT
  character (len=80) :: gftm_error_message='null'
  integer :: gftm_error_status=0
  
  ! Set threshold for gftm-NN execution versus full gftm calculation
  real    :: gftm_nn_max_error_in = 0.0

 ! CHECK <<<---------------------------------------------<<<
 !write(*,*) 'gftm_interface --> gftm_nn_max_error_in: ', gftm_nn_max_error_in

 
  
contains

  ! Dump LOCAL INTERFACE variables
  subroutine gftm_dump_local()

    implicit none

    integer :: ierr, i

    open(unit=1,file=trim(gftm_path_in)//trim(file_dump_local),&
         status='replace',iostat=ierr)

    if (gftm_geometry_flag_in == 2 ) then
       write(1,*) 'DUMP NOT SUPPORTED FOR FOURIER GEOMETRY'
       close(1)
       return
    endif

    write(1,*) '# input.gftm generated by gftm_dump_local()'
    write(1,*) '#'
    write(1,*) '# See https://fusion.gat.com/theory/gftminput'
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Control parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,50) 'UNITS', gftm_units_in
    write(1,20) 'NS',gftm_ns_in
    write(1,10) 'USE_TRANSPORT_MODEL',gftm_use_transport_model_in
    write(1,20) 'GEOMETRY_FLAG',gftm_geometry_flag_in
    write(1,10) 'USE_BPER',gftm_use_bper_in
    write(1,10) 'USE_BPAR',gftm_use_bpar_in
    write(1,10) 'USE_MHD_RULE',gftm_use_mhd_rule_in
    write(1,10) 'USE_BISECTION',gftm_use_bisection_in
    write(1,10) 'USE_INBOARD_DETRAPPED',gftm_use_inboard_detrapped_in
    write(1,10) 'USE_AVE_ION_GRID',gftm_use_ave_ion_grid_in
    write(1,20) 'SAT_RULE',gftm_sat_rule_in
    write(1,20) 'KYGRID_MODEL',gftm_kygrid_model_in
    write(1,20) 'XNU_MODEL',gftm_xnu_model_in
    write(1,20) 'VPAR_MODEL',gftm_vpar_model_in
    write(1,20) 'VPAR_SHEAR_MODEL',gftm_vpar_shear_model_in
    write(1,30) 'SIGN_BT',gftm_sign_bt_in
    write(1,30) 'SIGN_IT',gftm_sign_it_in
    write(1,30) 'KY',gftm_ky_in
    write(1,10) 'NEW_EIKONAL',gftm_new_eikonal_in
    write(1,30) 'VEXB',gftm_vexb_in
    write(1,30) 'VEXB_SHEAR',gftm_vexb_shear_in
    write(1,30) 'BETAE',gftm_betae_in
    write(1,30) 'XNUE',gftm_xnue_in
    write(1,30) 'ZEFF',gftm_zeff_in
    write(1,30) 'DEBYE',gftm_debye_in
    write(1,10) 'IFLUX',gftm_iflux_in
    write(1,20) 'IBRANCH',gftm_ibranch_in
    write(1,20) 'NMODES',gftm_nmodes_in
    write(1,20) 'NBASIS_MAX',gftm_nbasis_max_in
    write(1,20) 'NBASIS_MIN',gftm_nbasis_min_in
    write(1,20) 'NXGRID',gftm_nxgrid_in
    write(1,20) 'NKY',gftm_nky_in
    write(1,10) 'ADIABATIC_ELEC',gftm_adiabatic_elec_in
    write(1,30) 'ALPHA_MACH',gftm_alpha_mach_in
    write(1,30) 'ALPHA_E',gftm_alpha_e_in
    write(1,30) 'ALPHA_P',gftm_alpha_p_in
    write(1,30) 'ALPHA_QUENCH',gftm_alpha_quench_in
    write(1,30) 'ALPHA_ZF',gftm_alpha_zf_in
    write(1,30) 'XNU_FACTOR',gftm_xnu_factor_in
    write(1,30) 'DEBYE_FACTOR',gftm_debye_factor_in
    write(1,30) 'ETG_FACTOR',gftm_etg_factor_in
    write(1,30) 'RLNP_CUTOFF',gftm_rlnp_cutoff_in
    write(1,20) 'WRITE_WAVEFUNCTION_FLAG',gftm_write_wavefunction_flag_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Species vectors:'
    write(1,*) '#---------------------------------------------------'
    do i = 1,gftm_ns_in
        write(1,40) 'ZS_',i,gftm_zs_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'MASS_',i,gftm_mass_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'RLNS_',i,gftm_rlns_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'RLTS_',i,gftm_rlts_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'TAUS_',i,gftm_taus_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'AS_',i,gftm_as_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'VPAR_',i,gftm_vpar_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'VPAR_SHEAR_',i,gftm_vpar_shear_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'VNS_SHEAR_',i,gftm_vns_shear_in(i)
    enddo
    do i = 1,gftm_ns_in
        write(1,40) 'VTS_SHEAR_',i,gftm_vts_shear_in(i)
    enddo
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Gaussian width parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'WIDTH',gftm_width_in
    write(1,30) 'WIDTH_MIN',gftm_width_min_in
    write(1,10) 'FIND_WIDTH',gftm_find_width_in
    write(1,20) 'NWIDTH',gftm_nwidth_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Miller geometry parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'RMIN_LOC',gftm_rmin_loc_in
    write(1,30) 'RMAJ_LOC',gftm_rmaj_loc_in
    write(1,30) 'ZMAJ_LOC',gftm_zmaj_loc_in
    write(1,30) 'DRMINDX_LOC',gftm_drmindx_loc_in
    write(1,30) 'DRMAJDX_LOC',gftm_drmajdx_loc_in
    write(1,30) 'DZMAJDX_LOC',gftm_dzmajdx_loc_in
    write(1,30) 'Q_LOC',gftm_q_loc_in
    write(1,30) 'KAPPA_LOC',gftm_kappa_loc_in
    write(1,30) 'S_KAPPA_LOC',gftm_s_kappa_loc_in
    write(1,30) 'DELTA_LOC',gftm_delta_loc_in
    write(1,30) 'S_DELTA_LOC',gftm_s_delta_loc_in
    write(1,30) 'ZETA_LOC',gftm_zeta_loc_in
    write(1,30) 'S_ZETA_LOC',gftm_s_zeta_loc_in
    write(1,30) 'P_PRIME_LOC',gftm_p_prime_loc_in
    write(1,30) 'Q_PRIME_LOC',gftm_q_prime_loc_in
    write(1,30) 'KX0_LOC',gftm_kx0_loc_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# s-alpha geometry parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'RMIN_SA',gftm_rmin_sa_in
    write(1,30) 'RMAJ_SA',gftm_rmaj_sa_in
    write(1,30) 'Q_SA',gftm_q_sa_in
    write(1,30) 'SHAT_SA',gftm_shat_sa_in
    write(1,30) 'ALPHA_SA',gftm_alpha_sa_in
    write(1,30) 'XWELL_SA',gftm_xwell_sa_in
    write(1,30) 'THETA0_SA',gftm_theta0_sa_in
    write(1,20) 'B_MODEL_SA',gftm_b_model_sa_in
    write(1,20) 'FT_MODEL_SA',gftm_ft_model_sa_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Expert parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'DAMP_PSI',gftm_damp_psi_in
    write(1,30) 'DAMP_SIG',gftm_damp_sig_in
    write(1,30) 'PARK',gftm_park_in
    write(1,30) 'GHAT',gftm_ghat_in
    write(1,30) 'GCHAT',gftm_gchat_in
    write(1,30) 'WD_ZERO',gftm_wd_zero_in
    write(1,30) 'LINSKER_FACTOR',gftm_linsker_factor_in
    write(1,30) 'GRADB_FACTOR',gftm_gradB_factor_in
    write(1,30) 'FILTER',gftm_filter_in
    write(1,30) 'THETA_TRAPPED',gftm_theta_trapped_in
    write(1,30) 'NN_MAX_ERROR', gftm_nn_max_error_in

    close(1)

10  format(t2,a,'=',l1)
20  format(t2,a,'=',i3)
30  format(t2,a,'=',1pe12.5)
40  format(t2,a,i1,'=',1pe12.5)
50  format(t2,a,'=',a8)

  end subroutine gftm_dump_local

  ! Dump GLOBAL MODULE variables
  subroutine gftm_dump_global()

    use gftm_global

    implicit none

    integer :: ierr, i

    open(unit=1,file=trim(gftm_path_in)//'out.gftm.globaldump',&
         status='replace',iostat=ierr)
    if (gftm_geometry_flag_in == 2 ) then
       write(1,*) 'DUMP NOT SUPPORTED FOR FOURIER GEOMETRY'
       close(1)
       return
    endif

    write(1,*) '# input.gftm generated by gftm_dump_global()'
    write(1,*) '#'
    write(1,*) '# See https://fusion.gat.com/theory/gftminput'
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Control parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,50) 'UNITS',gftm_units_in
    write(1,20) 'NS',ns_in
    write(1,10) 'USE_TRANSPORT_MODEL',.TRUE.
    write(1,20) 'GEOMETRY_FLAG',igeo
    write(1,10) 'USE_BPER',use_bper_in
    write(1,10) 'USE_BPAR',use_bpar_in
    write(1,10) 'USE_MHD_RULE',use_mhd_rule_in
    write(1,10) 'USE_BISECTION',use_bisection_in
    write(1,10) 'USE_INBOARD_DETRAPPED',use_inboard_detrapped_in
    write(1,10) 'USE_AVE_ION_GRID',use_ave_ion_grid_in
    write(1,20) 'SAT_RULE',sat_rule_in
    write(1,20) 'KYGRID_MODEL',kygrid_model_in
    write(1,20) 'XNU_MODEL',xnu_model_in
    write(1,20) 'VPAR_MODEL',vpar_model_in
    write(1,20) 'VPAR_SHEAR_MODEL',vpar_shear_model_in 
    write(1,30) 'SIGN_BT',sign_bt_in
    write(1,30) 'SIGN_IT',sign_it_in
    write(1,30) 'KY',ky_in
    write(1,10) 'NEW_EIKONAL',new_eikonal_in
    write(1,30) 'VEXB',vexb_in
    write(1,30) 'VEXB_SHEAR',vexb_shear_in
    write(1,30) 'BETAE',betae_in
    write(1,30) 'XNUE',xnue_in
    write(1,30) 'ZEFF',zeff_in
    write(1,30) 'DEBYE',debye_in
    write(1,10) 'IFLUX',iflux_in
    write(1,20) 'IBRANCH',ibranch_in
    write(1,20) 'NMODES',nmodes_in
    write(1,20) 'NBASIS_MAX',nbasis_max_in
    write(1,20) 'NBASIS_MIN',nbasis_min_in
    write(1,20) 'NXGRID',nxgrid_in
    write(1,20) 'NKY',nky_in
    write(1,10) 'ADIABATIC_ELEC',adiabatic_elec_in
    write(1,30) 'ALPHA_MACH',alpha_mach_in
    write(1,30) 'ALPHA_E',alpha_e_in
    write(1,30) 'ALPHA_P',alpha_p_in
    write(1,30) 'ALPHA_QUENCH',alpha_quench_in
    write(1,30) 'ALPHA_ZF',alpha_zf_in
    write(1,30) 'XNU_FACTOR',xnu_factor_in
    write(1,30) 'DEBYE_FACTOR',debye_factor_in
    write(1,30) 'ETG_FACTOR',etg_factor_in
    write(1,30) 'RLNP_CUTOFF',rlnp_cutoff_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Species vectors:'
    write(1,*) '#---------------------------------------------------'
    do i = 1,ns_in
        write(1,40) 'ZS_',i,zs_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'MASS_',i,mass_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'RLNS_',i,rlns_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'RLTS_',i,rlts_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'TAUS_',i,taus_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'AS_',i,as_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'VPAR_',i,vpar_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'VPAR_SHEAR_',i,vpar_shear_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'VNS_SHEAR_',i,vns_shear_in(i)
    enddo
    do i = 1,ns_in
        write(1,40) 'VTS_SHEAR_',i,vts_shear_in(i)
    enddo
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Gaussian width parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'WIDTH',width_in
    write(1,30) 'WIDTH_MIN',width_min_in
    write(1,10) 'FIND_WIDTH',find_width_in
    write(1,20) 'NWIDTH',nwidth_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Miller geometry parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'RMIN_LOC',rmin_loc
    write(1,30) 'RMAJ_LOC',rmaj_loc
    write(1,30) 'ZMAJ_LOC',zmaj_loc
    write(1,30) 'DRMINDX_LOC',drmindx_loc
    write(1,30) 'DRMAJDX_LOC',drmajdx_loc
    write(1,30) 'DZMAJDX_LOC',dzmajdx_loc
    write(1,30) 'Q_LOC',q_loc
    write(1,30) 'KAPPA_LOC',kappa_loc
    write(1,30) 'S_KAPPA_LOC',s_kappa_loc
    write(1,30) 'DELTA_LOC',delta_loc
    write(1,30) 'S_DELTA_LOC',s_delta_loc
    write(1,30) 'ZETA_LOC',zeta_loc
    write(1,30) 'S_ZETA_LOC',s_zeta_loc
    write(1,30) 'P_PRIME_LOC',p_prime_loc
    write(1,30) 'Q_PRIME_LOC',q_prime_loc
    write(1,30) 'KX0_LOC',kx0_loc
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# s-alpha geometry parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'RMIN_SA',rmin_sa
    write(1,30) 'RMAJ_SA',rmaj_sa
    write(1,30) 'Q_SA',q_sa
    write(1,30) 'SHAT_SA',shat_sa
    write(1,30) 'ALPHA_SA',alpha_sa
    write(1,30) 'XWELL_SA',xwell_sa
    write(1,30) 'THETA0_SA',theta0_sa
    write(1,20) 'B_MODEL_SA',b_model_sa
    write(1,20) 'FT_MODEL_SA',ft_model_sa
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Expert parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'DAMP_PSI',damp_psi_in
    write(1,30) 'DAMP_SIG',damp_sig_in
    write(1,30) 'PARK',park_in
    write(1,30) 'GHAT',ghat_in
    write(1,30) 'GCHAT',gchat_in
    write(1,30) 'WD_ZERO',wd_zero_in
    write(1,30) 'LINSKER_FACTOR',linsker_factor_in
    write(1,30) 'GRADB_FACTOR',gradB_factor_in
    write(1,30) 'FILTER',filter_in
    write(1,30) 'THETA_TRAPPED',theta_trapped_in
    write(1,30) 'NN_MAX_ERROR', nn_max_error_in

    close(1)

10  format(t2,a,'=',l1)
20  format(t2,a,'=',i3)
30  format(t2,a,'=',1pe12.5)
40  format(t2,a,i1,'=',1pe12.5)
50  format(t2,a,'=',a8)

  end subroutine gftm_dump_global

end module gftm_interface
