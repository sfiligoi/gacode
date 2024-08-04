      MODULE gftm_max_dimensions
!
! global dimensions for shared arrays
!
      IMPLICIT NONE
      SAVE
!
      INTEGER, PARAMETER :: nb=32
      INTEGER, PARAMETER :: nxm=2*(nb+1)-1
      INTEGER, PARAMETER :: nsm=12, nt0=40
      INTEGER, PARAMETER :: nkym=512
      INTEGER, PARAMETER :: maxmodes=16
      INTEGER, PARAMETER :: max_ELITE=700
      INTEGER, PARAMETER :: max_fourier = 24
      INTEGER, PARAMETER :: ms = 128  ! ms needs to be divisible by 8
      INTEGER, PARAMETER :: max_plot =18*ms/8+1
      INTEGER, PARAMETER :: num = 13
      INTEGER, PARAMETER :: nem = 7
!
      END MODULE gftm_max_dimensions
!
      MODULE gftm_dimensions
      USE gftm_max_dimensions
!
! dimensions determined by inputs
! 
      IMPLICIT NONE
      SAVE
!
! global dimensions
!
      INTEGER nx,nbasis,nbasis_max,ns0,ns,nky
      INTEGER :: nu=3, ne=2, nz, nune, nphase, ntot
!
      END MODULE gftm_dimensions
!
      MODULE gftm_global
      USE gftm_max_dimensions
!
!  global controls for the gftm driver routine
!
      IMPLICIT NONE
      SAVE
! global constants
      REAL :: pi
      REAL :: pi_2
      REAL :: sqrt_pi
      REAL :: sqrt_two
      COMPLEX :: xi=(0.0,1.0)
! internal flow control switches
      LOGICAL :: new_eikonal_in=.TRUE.
      LOGICAL :: new_start=.TRUE.
      LOGICAL :: new_matrix=.TRUE.
      LOGICAL :: new_geometry=.TRUE.
      LOGICAL :: new_width=.TRUE.
      LOGICAL :: new_kyspectrum=.TRUE.
      LOGICAL :: gauher_uncalled=.TRUE.
      LOGICAL :: gauss_hermite_uncalled=.TRUE.
      LOGICAL :: eikonal_unsaved=.TRUE.
      INTEGER :: igeo=1
      LOGICAL :: use_default_species=.TRUE.
      INTEGER,DIMENSION(8) :: trace_path=0
      LOGICAL,EXTERNAL :: gftm_isnan
      LOGICAL,EXTERNAL :: gftm_isinf
! input units switch
      CHARACTER (len=8) :: units_in = 'GYRO'
! Input Gaussian width
      REAL :: width_in=1.65
      REAL :: width_min_in=0.3
      INTEGER :: nwidth_in=21
      LOGICAL :: find_width_in=.TRUE.
! Input kys
      REAL :: ky_in=0.3
! Input species 
      INTEGER :: ns_in=2, nstotal_in = 2
      REAL,DIMENSION(nsm) :: mass_in
      REAL,DIMENSION(nsm) :: zs_in
! input switches
      LOGICAL :: iflux_in=.TRUE.
      LOGICAL :: use_bper_in=.FALSE.
      LOGICAL :: use_bpar_in=.FALSE.
      LOGICAL :: use_mhd_rule_in=.TRUE.
      LOGICAL :: use_bisection_in=.TRUE.
      LOGICAL :: use_inboard_detrapped_in=.FALSE.
      LOGICAL :: use_ave_ion_grid_in=.FALSE.
      INTEGER :: ibranch_in=-1
      INTEGER :: nmodes_in=2
      INTEGER :: nbasis_max_in=4
      INTEGER :: nbasis_min_in=2
      INTEGER :: nxgrid_in=16
      INTEGER :: nky_in=12
      INTEGER :: mainion=2
! input rare switches
      REAL :: theta_trapped_in=0.7
      REAL :: wdia_trapped_in=0.0
      REAL :: park_in=1.0
      REAL :: ghat_in=1.0
      REAL :: gchat_in=1.0
      REAL :: wd_zero_in=0.1
      REAL :: Linsker_factor_in=0.0
      REAL :: gradB_factor_in=0.0
      REAL :: filter_in=2.0
      REAL :: damp_psi_in = 0.0
      REAL :: damp_sig_in = 0.0
! Input model paramaters
      LOGICAL :: adiabatic_elec_in=.FALSE.
      REAL :: alpha_e_in=1.0
      REAL :: alpha_p_in=1.0
      REAL :: alpha_mach_in=0.0
      REAL :: alpha_quench_in=0.0
      REAL :: alpha_zf_in = 1.0
      REAL :: xnu_factor_in=1.0
      REAL :: debye_factor_in=1.0
      REAL :: etg_factor_in=1.25
      REAL :: rlnp_cutoff_in = 18.0
      INTEGER :: sat_rule_in=0
      INTEGER :: kygrid_model_in=1
      INTEGER :: xnu_model_in=2
      INTEGER :: vpar_model_in=0
      INTEGER :: vpar_shear_model_in=1    
      REAL :: alpha_n_in =0.0  !not used
      REAL :: alpha_t_in =0.0  !not used
! Input signs
      REAL :: sign_Bt_in = 1.0
      REAL :: sign_It_in = 1.0
! Input field gradients
      REAL,DIMENSION(nsm) :: rlns_in=0.0
      REAL,DIMENSION(nsm) :: rlts_in=0.0
      REAL,DIMENSION(nsm) :: vpar_shear_in=0.0
      REAL :: vexb_shear_in=0.0
! Input profile shear
      REAL,DIMENSION(nsm) :: vns_shear_in=0.0
      REAL,DIMENSION(nsm) :: vts_shear_in=0.0
! Input field averages
      REAL,DIMENSION(nsm) :: as_in=0.0
      REAL,DIMENSION(nsm) :: taus_in=0.0
      REAL,DIMENSION(nsm) :: vpar_in=0.0
      REAL :: vexb_in = 0.0
      REAL :: betae_in=0.0
      REAL :: xnue_in=0.0
      REAL :: zeff_in=1.0
      REAL :: debye_in=0.0
! Shifted circle (s-alpha) inputs
      REAL :: rmin_sa=0.5
      REAL :: rmaj_sa=3.0
      REAL :: q_sa=2.0
      REAL :: shat_sa=1.0
      REAL :: alpha_sa=0.0
      REAL :: xwell_sa=0.0
      REAL :: theta0_sa=0.0
! Shifted circle flags
      INTEGER :: b_model_sa=1
      INTEGER :: ft_model_sa=1
! Miller inputs
      REAL :: rmin_loc=0.5
      REAL :: rmaj_loc=3.0
      REAL :: zmaj_loc=0.0
      REAL :: q_loc=2.0
!      REAL :: shat_loc=1.0
      REAL :: dlnpdr_loc=0.0
      REAL :: drmindx_loc=1.0
      REAL :: drmajdx_loc=0.0
      REAL :: dzmajdx_loc=0.0
      REAL :: kappa_loc=1.0
      REAL :: s_kappa_loc=0.0
      REAL :: delta_loc=0.0
      REAL :: s_delta_loc=0.0
      REAL :: zeta_loc=0.0
      REAL :: s_zeta_loc=0.0
      REAL :: p_prime_loc=0.0 
      REAL :: q_prime_loc=16.0
      REAL :: kx0_loc = 0.0
! Fourier geometry inputs
      INTEGER :: nfourier_in        = 16
      REAL :: q_fourier_in          = 2.0
      REAL :: q_prime_fourier_in    = 16.0
      REAL :: p_prime_fourier_in    = 0.0
      REAL :: fourier_in(8,0:max_fourier) = 0.0
! ELITE geometry inputs
      INTEGER :: n_ELITE
      REAL :: R_ELITE(0:max_ELITE)
      REAL :: Z_ELITE(0:max_ELITE)
      REAL :: Bp_ELITE(0:max_ELITE)
      REAL :: q_ELITE
      REAL :: q_prime_ELITE
      REAL :: p_prime_ELITE
! global variables
      REAL :: ft=0.5
      REAL :: ft_min=0.01
      REAL :: ft_test
      REAL :: modB_test
      REAL :: modB_min
      REAL :: ky=0.3
      REAL :: R_unit=3.0
      REAL :: q_unit=2.0
      REAL :: B_unit=1.0
      REAL :: Rmaj_input = 3.0
      REAL :: q_in = 2.0
      REAL :: rmin_input=0.5
      REAL,DIMENSION(maxmodes) :: gamma_reference_kx0 =0.0
      REAL,DIMENSION(maxmodes) :: freq_reference_kx0 =0.0
      REAL,DIMENSION(2,nkym,maxmodes) :: eigenvalue_first_pass =0.0
      REAL :: pol=1.0
      REAL :: U0=0.0
      REAL :: kx0=0.0
      REAL :: kx0_e=0.0
      REAL :: kx0_p = 0.0
      REAL :: midplane_shear=1.0
      REAL :: kx0_factor=1.0
      REAL :: rho_ion=1.0
      REAL :: rho_e=1.0
! output
      COMPLEX,DIMENSION(3,nb) :: field_weight_QL_out=0.0
      COMPLEX,DIMENSION(maxmodes,3,nb) :: field_weight_out=0.0
      COMPLEX,DIMENSION(maxmodes,3,max_plot) :: plot_field_out=0.0
      REAL,DIMENSION(max_plot) :: plot_angle_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: particle_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: energy_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: stress_par_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: stress_tor_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: exchange_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: N_QL_out=0.0,T_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: U_QL_out=0.0,Q_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: n_bar_out=0.0,t_bar_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: u_bar_out=0.0,q_bar_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: Ns_Ts_phase_out=0.0
      REAL,DIMENSION(nsm) :: particle_flux_out=0.0,energy_flux_out=0.0
      REAL,DIMENSION(nsm) :: exchange_out=0.0
      REAL,DIMENSION(nsm) :: stress_par_out=0.0,stress_tor_out=0.0
      REAL,DIMENSION(maxmodes) :: gamma_out=0.0,freq_out=0.0
      REAL,DIMENSION(maxmodes) :: v_QL_out=0.0,a_par_QL_out=0.0,b_par_QL_out=0.0
      REAL,DIMENSION(maxmodes) :: phi_bar_out=0.0,v_bar_out=0.0
      REAL,DIMENSION(maxmodes) :: a_par_bar_out=0.0,b_par_bar_out=0.0
      REAL,DIMENSION(maxmodes) :: wd_bar_out=0.0,b0_bar_out=0.0
      REAL,DIMENSION(maxmodes) :: ne_te_phase_out=0.0
      REAL,DIMENSION(maxmodes) :: kx_bar_out=0.0,kpar_bar_out=0.0
      REAL,DIMENSION(maxmodes) :: modB_bar_out=0.0
      REAL,DIMENSION(nsm) :: n_bar_sum_out=0.0,t_bar_sum_out=0.0
      REAL,DIMENSION(nsm) :: q_low_out=0.0
      REAL,DIMENSION(4,nkym,maxmodes) :: field_spectrum_out=0.0
      REAL,DIMENSION(4,nkym,maxmodes) :: QL_field_spectrum_out=0.0
      REAL,DIMENSION(4,nsm,nkym,maxmodes) :: intensity_spectrum_out=0.0
      REAL,DIMENSION(4,nsm,nkym,maxmodes) :: QL_intensity_spectrum_out=0.0
      REAL,DIMENSION(5,nsm,nkym,maxmodes) :: flux_spectrum_out=0.0
      REAL,DIMENSION(5,nsm,nkym,maxmodes) :: QL_flux_spectrum_out=0.0
      REAL,DIMENSION(2,nkym,maxmodes) :: eigenvalue_spectrum_out=0.0
      REAl,DIMENSION(nkym,maxmodes) :: ne_te_phase_spectrum_out=0.0
      REAl,DIMENSION(nsm,nkym,maxmodes) :: nsts_phase_spectrum_out=0.0
      REAL,DIMENSION(nkym) :: spectral_shift_out=0.0
      REAL,DIMENSION(nkym) :: ave_p0_spectrum_out=0.0
      REAL,DIMENSION(nkym) :: width_out=0.0
      REAL,DIMENSION(nsm,5,5) :: diff_out=0.0
      REAL,DIMENSION(nsm,5) :: conv_out=0.0
      REAL,DIMENSION(nsm,5) :: flux_out=0.0
      REAL,DIMENSION(nsm) :: delta_out
      REAL :: Vzf_out = 0.0
      REAL :: kymax_out = 0.0
      REAL :: phi_bar_sum_out=0.0
      REAL :: v_bar_sum_out=0.0
      REAL :: gamma_nb_min_out=0.0
      REAL :: B2_ave_out=1.0
      REAL :: R2_ave_out=1.0
      REAL :: B_ave_out=1.0
      REAL :: Bt_ave_out=1.0
      REAL :: Bp0_out = 1.0
      REAL :: RBt_ave_out=1.0
      REAL :: Grad_r_ave_out=1.0
      REAL :: grad_r0_out=1.0
      REAL :: B_geo0_out = 1.0
      REAL :: Bt0_out = 1.0
      REAL :: SAT_geo0_out=1.0
      REAL :: SAT_geo1_out=1.0
      REAL :: SAT_geo2_out=1.0
      REAL :: kx_geo0_out=1.0
      REAL :: DM_out = 0.25
      REAL :: DR_out = 0.0
      REAL :: Bref_out = 1.0
      REAL :: ave_p0_out = 1.0
      INTEGER :: nmodes_out = 2
      INTEGER :: nfields_out = 1
      INTEGER :: jmax_out = 0
      character (len=80) :: error_msg='null' 
! NN activation parameters (thresholds)  
      REAL :: nn_max_error_in = -1.0
      LOGICAL :: valid_nn = .FALSE.
    
     
!
      END MODULE gftm_global
!-------------------------------------------------
      MODULE gftm_hermite
! hermite basis functions and x-grid
      USE gftm_max_dimensions
      IMPLICIT NONE
      SAVE
!
      REAL :: x(nxm),wx(nxm),h(nxm,nxm)
!
      END MODULE gftm_hermite
!
      MODULE gftm_species
! species parameters
!      USE gftm_max_dimensions
      IMPLICIT NONE
!
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ei_exch, resist
      REAL,ALLOCATABLE,DIMENSION(:) :: zs, mass, vs, fts
      REAL,ALLOCATABLE,DIMENSION(:) :: rlts, rlns
      REAL,ALLOCATABLE,DIMENSION(:) :: as, taus
      REAL,ALLOCATABLE,DIMENSION(:) :: vpar_s, vpar_shear_s
      REAL :: xnue_s, vexb_shear_s, ky_s
!
      END MODULE gftm_species
!-------------------------------------------------
      MODULE gftm_kyspectrum
!
! ky spectrum for computing total fluxes and intensities
!
      USE gftm_dimensions
      IMPLICIT NONE
      SAVE

!
      REAL,DIMENSION(nkym) :: ky_spectrum
      REAL,DIMENSION(nkym) :: dky_spectrum
!
      END MODULE gftm_kyspectrum
!
!-------------------------------------------------
      MODULE gftm_eigen
!
! eigenvalues, eigenvectors and fluxes
!
!      USE gftm_dimensions
      IMPLICIT NONE
!
      INTEGER matz
      INTEGER,ALLOCATABLE,DIMENSION(:) :: grow_index
      REAL,ALLOCATABLE,DIMENSION(:) :: fv1, fv2, fv3
      REAL,ALLOCATABLE,DIMENSION(:) :: rr, ri
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ar, ai
      REAL,ALLOCATABLE,DIMENSION(:,:) :: vr, vi
      COMPLEX :: we
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: amat, bmat
      COMPLEX,ALLOCATABLE,DIMENSION(:) :: he, alpha, beta
      COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: hetot
!
      END MODULE gftm_eigen
!
!-------------------------------------------------
!
      MODULE gftm_weight
!  
!  quasilinear weights and eigenfunction averages
!
      USE gftm_dimensions
      IMPLICIT NONE
!
      REAL :: N_weight(nsm),T_weight(nsm)
      REAl :: U_weight(nsm),Q_weight(nsm)
      REAL :: phi_weight,a_par_weight,b_par_weight,v_weight
      REAL :: Ne_Te_phase,Ne_Te_cos,Ne_Te_sin
      REAL :: Ns_Ts_phase(nsm),Ns_Ts_cos,Ns_Ts_sin
      REAL :: wd_bar,b0_bar,modB_bar,kx_bar,kpar_bar
!      
      END MODULE gftm_weight
!
!-------------------------------------------------
!      
     MODULE gftm_xgrid
!
! functions on the x-grid
!
      USE gftm_max_dimensions
      IMPLICIT NONE
      SAVE
!
      REAL,DIMENSION(nxm) :: wdx, wdpx, b0x, b2x, kxx
      REAL,DIMENSION(nxm) :: cx_tor_par, cx_tor_per
      REAL,DIMENSION(nxm) :: cx_par_par
      REAL,DIMENSION(nxm) :: p0x, Bx
      INTEGER,DIMENSION(nkym) :: mask_save
      REAL,DIMENSION(nkym) :: gamma_nb_min_save
      REAL,DIMENSION(nkym) :: width_save, ft_save
      REAL,DIMENSION(nkym) :: R_unit_save, q_unit_save
      REAL,DIMENSION(nkym,nxm) :: wdx_save, b0x_save
      REAL,DIMENSION(nkym,nxm) :: cx_par_par_save
      REAL,DIMENSION(nkym,nxm) :: cx_tor_par_save
      REAL,DIMENSION(nkym,nxm) :: cx_tor_per_save
      REAL,DIMENSION(nkym,nxm) :: b2x_save, kxx_save
! gftm2 additions
      REAL,DIMENSION(nsm,nem,nxm) :: pe1j0sx,pe2j0sx,pe1j1sx,pe2j1sx
!
      END MODULE gftm_xgrid
!
!-------------------------------------------------
!
      MODULE gftm_sgrid
!
! functions on the s-grid
!
      USE gftm_max_dimensions
      IMPLICIT NONE
      SAVE
!
!---------------------------------------------------------------
! s_grid.m
! 
! PURPOSE:
!  include file which contains the data for mercier_luc.f :
! INPUT
!  R, Z and Bp on the s-grid with ms equally space intervals of length ds
!  as defined by the polidal co-ordinate in the mercier-luc system. 
!  B_unit = B0 (r/rho) dr/drho where B0 = magnetic field on axis
!  Rmaj_s = major radius at magnetic axis
!  rmin_s = minor radius of flux surface 
!  q_s = local flux surface safety factor 
!  q_prime_s = dq/dpsi
!  p_prime_s = dp/dpsi 
!
! OUTPUT
!  costheta_geo, sintheta_geo, costheta_p_geo, pk_geo, epsl_geo, qrat_geo
!  kyoky_geo, b_geo
!
! 24 June 05: gms
!  this version for gftm separated miller and mercier-luc components
!  in GKS and GYRO versions
!---------------------------------------------------------------
! INPUT
      REAL,DIMENSION(0:ms) :: R, Z, Bp
      REAL :: ds, Ls
      REAL :: Rmaj_s, Zmaj_s,rmin_s, q_s
      REAL :: p_prime_s, q_prime_s
      REAL :: p_prime_zero_s
      REAL :: betae_s, debye_s
! OUTPUT
      REAL,DIMENSION(0:ms) :: costheta_geo, sintheta_geo
      REAL,DIMENSION(0:ms) :: costheta_p_geo, s_p
      REAL,DIMENSION(0:ms) :: pk_geo, epsl_geo, qrat_geo
      REAL,DIMENSION(0:ms) :: kxoky_geo, b_geo, t_s
      REAL,DIMENSION(0:ms) :: S_prime, kx_factor, y
      REAL :: f,ff_prime
!
      END MODULE gftm_sgrid
!
!-------------------------------------------------
!
      MODULE gftm_coeff

 !  ave_theta
      REAL :: gradB
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_kx
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_c_tor_par, ave_c_tor_per
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_c_par_par
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_wdpar
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_wdper
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_gradB, ave_lnB
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_b0, ave_b0inv
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_kpar, ave_modkpar
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_p0, ave_p0inv
      REAL,ALLOCATABLE,DIMENSION(:,:) :: ave_bp, ave_bpinv
!
      END MODULE gftm_coeff
!
!-----------------------------------------------------
!
       MODULE gftm_GFS
       REAL,ALLOCATABLE,DIMENSION(:,:) :: matmirror
       REAL,ALLOCATABLE,DIMENSION(:,:) :: matu, matuc,matdu, matuu
       REAL,ALLOCATABLE,DIMENSION(:,:) :: mate, matde
       REAL,ALLOCATABLE,DIMENSION(:,:) :: phib, psib, sigb
       REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: pe1j0phib, pe1j0psib, pe1j1sigb
       REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: pe2j0phib, pe2j0psib, pe2j1sigb
       REAL,ALLOCATABLE,DIMENSION(:,:) :: pu1pe1PSI, pu1pe2PSI
       REAL,ALLOCATABLE,DIMENSION(:,:) :: pu3pe1PSI,pu2pe1PSI
       REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: p0invpe1j0s, b0invpe1j0s
       COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: mateq, mats, matb
       COMPLEX,ALLOCATABLE,DIMENSION(:) :: zomega
!
       INTEGER,ALLOCATABLE,DIMENSION(:) :: di,de
       INTEGER,ALLOCATABLE,DIMENSION(:) :: ipiv,ipive
       COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: emat,zmat, gmat
       COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: gmatinv
!
       COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: pe1j0phi, pe1j0psi, pe1j1sig
       COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: pe2j0phi, pe2j0psi, pe2j1sig
       COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: pwns,pwps,pwtpars,pwtpers,pwes
       COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: hwns,hwps,hwtpars,hwtpers,hwes,hes
       COMPLEX,ALLOCATABLE,DIMENSION(:) :: phi,psi,sig
       COMPLEX,ALLOCATABLE,DIMENSION(:,:)  :: PSIpu1pe1, PSIpu2pe1
       COMPLEX,ALLOCATABLE,DIMENSION(:,:)  :: PSIpu3pe1, PSIpu1pe2
!
       COMPLEX,ALLOCATABLE,DIMENSION(:,:) :: rwn,rwp,rwtpar,rwtper,rwe,fluxe
!
       END MODULE gftm_GFS
!
!-------------------------------------------------
!

!-----------------------------------------------------
!
      MODULE gftm_gyro_average
!
!     integrals of gyro average operators
!
      USE gftm_max_dimensions
      IMPLICIT NONE
!
      REAL,DIMENSION(nem) :: pe1j0=0.0
      REAL,DIMENSION(nem) :: pe2j0=0.0
      REAL,DIMENSION(nem) :: pe1j1=0.0
      REAL,DIMENSION(nem) :: pe2j1=0.0
      REAL,DIMENSION(nsm,nem,nb,nb) :: mat_pe1j0s=0.0
      REAL,DIMENSION(nsm,nem,nb,nb) :: mat_pe2j0s=0.0
      REAL,DIMENSION(nsm,nem,nb,nb) :: mat_pe1j1s=0.0
      REAL,DIMENSION(nsm,nem,nb,nb) :: mat_pe2j1s=0.0
 !
      END MODULE gftm_gyro_average
!-----------------------------------------------------
!
      MODULE gftm_velocity_matrix
!
!     integrals of gyro average operators
!
      USE gftm_max_dimensions
      IMPLICIT NONE
!
      REAL,DIMENSION(num,num) :: mat_upar
      REAL,DIMENSION(num,num) :: mat_uparc
      REAL,DIMENSION(num,num) :: mat_dupar
!
      REAL,DIMENSION(nem,nem) :: mat_eper
      REAL,DIMENSION(nem,nem) :: mat_deper
 !
      END MODULE gftm_velocity_matrix
!-----------------------------------------------------
!
