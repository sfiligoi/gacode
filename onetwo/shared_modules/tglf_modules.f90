      MODULE tglf_dimensions
!
! global dimensions for shared arrays
!
      IMPLICIT NONE
      SAVE
!
      INTEGER, PARAMETER:: nb=24,nxm=47
      INTEGER, PARAMETER:: nsm=3,nt0=20
      INTEGER, PARAMETER :: neq = 15*nsm,iar=neq*nb
      INTEGER, PARAMETER :: nkym=50
      INTEGER, PARAMETER :: maxmodes=4
! dimensions determined by inputs
      INTEGER nx,nbasis,ns0,ns
! global constants
      REAL :: pi
      REAL :: pi_2
      REAL :: sqrt_pi
      REAL :: sqrt_two
!
      END MODULE tglf_dimensions
 !
      MODULE tglf_interface
      USE tglf_dimensions
!
!  global controls for the tglf driver routine
!
      IMPLICIT NONE
      SAVE
! internal flow control switches
      LOGICAL :: new_eikonal_in=.TRUE.
      LOGICAL :: new_start=.TRUE.
      LOGICAL :: new_matrix=.TRUE.
      LOGICAL :: new_geometry=.TRUE.
      LOGICAL :: new_width=.TRUE.
      LOGICAL :: new_kyspectrum=.TRUE.
      LOGICAL :: species_defined=.FALSE.
      LOGICAL :: eikonal_unsaved=.TRUE.
      INTEGER :: igeo=0
! Input global variables
      LOGICAL :: find_width_in=.FALSE.
      INTEGER :: nwidth_in=21
      REAL :: width_in=1.65
      REAL :: width_min_in=0.15
      REAL :: ky_in=0.3
      INTEGER :: nky_in=15
      INTEGER :: nmodes_in=2
! Input species 
      INTEGER :: ns_in=2
      REAL :: mass_in(nsm)
      REAL :: zs_in(nsm)
! input switches
      LOGICAL :: iflux_in=.TRUE.
      LOGICAL :: use_bper_in=.FALSE.
      LOGICAL :: use_bpar_in=.FALSE.
      LOGICAL :: use_bisection_in=.FALSE.
      INTEGER :: ibranch_in=0
      INTEGER :: nbasis_max_in=4
      INTEGER :: nbasis_min_in=1
      INTEGER :: nxgrid_in=16
! input rare switches
      REAL :: park_in=1.0
      REAL :: ghat_in=1.0
      REAL :: gchat_in=1.0
      REAL :: wd_zero_in=0.1
      REAL :: Linsker_factor_in=0.0
      REAL :: gradB_factor_in=0.0
      REAL :: filter_in=2.0
      REAL :: x_psi_in=0.0
! Input model paramaters
      LOGICAL :: adiabatic_elec_in=.FALSE.
      REAL :: alpha_p_in=0.0
      REAL :: alpha_e_in=0.0
      REAL :: theta_trapped_in=0.7
      REAL :: xnu_factor_in=1.0
      REAL :: debye_factor_in=1.0
      REAL :: etg_factor_in=4.0
! Input field gradients
      REAL :: rlns_in(nsm)
      REAL :: rlts_in(nsm)
      REAL :: vexb_shear_in=0.0
      REAL :: vpar_shear_in=0.0
! Input field averages
      REAL :: as_in(nsm)
      REAL :: taus_in(nsm)
      REAL :: betae_in=0.0
      REAL :: xnuei_in=0.0
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
      REAL :: rmin_loc=3.0
      REAL :: rmaj_loc=0.5
      REAL :: q_loc=2.0
      REAL :: shat_loc=1.0
      REAL :: dlnpdr_loc
      REAL :: shift_loc=0.0
      REAL :: kappa_loc=1.0
      REAL :: s_kappa_loc=0.0
      REAL :: delta_loc=0.0
      REAL :: s_delta_loc=0.0
      REAL :: p_prime_loc=0.0 
      REAL :: q_prime_loc=4.0
! global variables
      REAL :: ft
      REAL :: ft_min=0.01
      REAL :: ky
      REAL :: R_unit
      REAL :: q_unit
! output
      REAL :: gamma_out(maxmodes),freq_out(maxmodes)
      REAL :: gamma_net_out(maxmodes)
      REAL :: wd_bar_out(maxmodes),phi_QL_out(maxmodes)
      REAl :: particle_QL_out(maxmodes,nsm,3)
      REAL :: energy_QL_out(maxmodes,nsm,3)
      REAL :: exchange_QL_out(maxmodes,nsm,3)
      REAL :: N_QL_out(maxmodes,nsm),T_QL_out(maxmodes,nsm)
      REAL :: phi_bar_out(maxmodes),v_bar_out(maxmodes)
      REAL :: n_bar_out(maxmodes,nsm),t_bar_out(maxmodes,nsm)
      REAL :: particle_flux_out(nsm,3),energy_flux_out(nsm,3)
      REAL :: exchange_out(nsm,3)
      REAL :: stress_tor_out(nsm,3),stress_par_out(nsm,3)
      REAL :: phi_bar_sum_out
      REAL :: v_bar_sum_out=0.0
      REAL :: n_bar_sum_out(nsm),t_bar_sum_out(nsm)
      REAL :: q_low_out(nsm)
      REAL :: gamma_nb_min_out
!
      END MODULE tglf_interface
!------------------------------------------------- 
      MODULE tglf_closure
!
!  tglf closure coefficients
!
      IMPLICIT NONE
      SAVE
! ft = 1 coefficients
      REAL :: v1_r,v2_r,v3_r,v4_r,v5_r
      REAL :: v1_i,v2_i,v3_i,v4_i,v5_i
      REAL :: v6_r,v7_r,v8_r,v9_r,v10_r
      REAL :: v6_i,v7_i,v8_i,v9_i,v10_i
      REAL :: vb1_r,vb2_r,vb3_r,vb4_r,vb5_r
      REAL :: vb1_i,vb2_i,vb3_i,vb4_i,vb5_i
      REAL :: vb6_r,vb7_r,vb8_r,vb9_r,vb10_r
      REAL :: vb6_i,vb7_i,vb8_i,vb9_i,vb10_i
! ft < 1 coefficients
      REAL :: u1_r,u2_r,u3_r,u4_r,u5_r
      REAL :: u1_i,u2_i,u3_i,u4_i,u5_i
      REAL :: u6_r,u7_r,u8_r,u9_r,u10_r
      REAL :: u6_i,u7_i,u8_i,u9_i,u10_i
      REAL :: ub1_r,ub2_r,ub3_r,ub4_r,ub5_r
      REAL :: ub1_i,ub2_i,ub3_i,ub4_i,ub5_i
      REAL :: ub6_r,ub7_r,ub8_r,ub9_r,ub10_r
      REAL :: ub6_i,ub7_i,ub8_i,ub9_i,ub10_i
! parallel coefficients
      REAL :: dpar_HP,dper_HP,bpar_HP,bper_HP
      REAL :: bper_DH,dper_DH
      COMPLEX :: dhr13,dgr13
      REAL :: d1,d3,d33
      REAL :: b1,b3,b33
!
      END MODULE tglf_closure
!-------------------------------------------------
      MODULE tglf_hermite
! hermite basis functions and x-grid
      USE tglf_dimensions
      IMPLICIT NONE
      SAVE
!
      REAL :: x(nxm),wx(nxm),h(nb,nxm)
!
      END MODULE tglf_hermite  
!
      MODULE tglf_species
! species parameters
      USE tglf_dimensions
      IMPLICIT NONE
      SAVE
!
      REAL :: ei_exch(nsm,nsm),resist(nsm,nsm)
      REAL :: zs(nsm),mass(nsm),taus(nsm),vs(nsm)
      REAL :: rlts(nsm),rlns(nsm),as(nsm)
! 
      END MODULE tglf_species  
!-------------------------------------------------
      MODULE tglf_kyspectrum
!
! ky spectrum for computing total fluxes and intensities
!
      USE tglf_dimensions
      IMPLICIT NONE
      SAVE
!
      INTEGER nky
!
      REAL :: ky_spectrum(nkym)
      REAL :: dky_spectrum(nkym)
!
      END MODULE tglf_kyspectrum
!
!-------------------------------------------------
      MODULE tglf_eigen
!
! eigenvalues, eigenvectors and fluxes
!
      USE tglf_dimensions
      IMPLICIT NONE
      SAVE
!
      INTEGER matz,iur,nroot
      REAL :: fv1(iar),fv2(iar),fv3(iar)
      REAL :: rr(iar), ri(iar)
      REAL :: ar(iar,iar), ai(iar,iar)
      REAL :: vr(iar,iar), vi(iar,iar)
      COMPLEX :: v(iar),eigenvalue
      COMPLEX :: amat(iar,iar),bmat(iar,iar)
      COMPLEX :: alpha(iar),beta(iar)
!
      END MODULE tglf_eigen
!
!-------------------------------------------------
!
      MODULE tglf_xgrid
!
! functions on the x-grid
!
      USE tglf_dimensions
      IMPLICIT NONE
      SAVE
!
      REAL :: hxn(nsm,nxm),hxp1(nsm,nxm),hxp3(nsm,nxm)
      REAL :: hxr11(nsm,nxm),hxr13(nsm,nxm),hxr33(nsm,nxm)
      REAL :: hxw113(nsm,nxm),hxw133(nsm,nxm),hxw333(nsm,nxm)
      REAL :: gxn(nsm,nxm),gxp1(nsm,nxm),gxp3(nsm,nxm)
      REAL :: gxr11(nsm,nxm),gxr13(nsm,nxm),gxr33(nsm,nxm)
      REAL :: gxw113(nsm,nxm),gxw133(nsm,nxm),gxw333(nsm,nxm)
      INTEGER :: mask_save(nkym)
      REAL :: wdx(nxm),b0x(nxm)
      REAL :: p0x(nxm),Bx(nxm)
      REAL :: gamma_nb_min_save(nkym)
      REAL :: width_save(nkym),ft_save(nkym)
      REAL :: R_unit_save(nkym),q_unit_save(nkym)
      REAL :: wdx_save(nkym,nxm),b0x_save(nkym,nxm)
!
      END MODULE tglf_xgrid
!
!-------------------------------------------------
!
      MODULE tglf_sgrid
!
! functions on the s-grid
!
      USE tglf_dimensions
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
!  this version for TGLF separated miller and mercier-luc components
!  in GKS and GYRO versions
!---------------------------------------------------------------
      INTEGER, PARAMETER :: ms = 128  ! ms needs to be even
! INPUT
      REAL :: R(0:ms), Z(0:ms), Bp(0:ms)
      REAL :: ds, Ls
      REAL :: Rmaj_s, rmin_s, q_s
      REAL :: p_prime_s, q_prime_s
      REAL :: p_prime_zero_s
! OUTPUT
      REAL :: costheta_geo(0:ms),sintheta_geo(0:ms)
      REAL :: costheta_p_geo(0:ms),s_p(0:ms)
      REAL :: pk_geo(0:ms),epsl_geo(0:ms),qrat_geo(0:ms)
      REAL :: kxoky_geo(0:ms), b_geo(0:ms),t_s(0:ms)
      REAL :: S_prime(0:ms),kx_factor(0:ms),y(0:ms)
      REAL :: f,ff_prime
!
      END MODULE tglf_sgrid  
!
!-------------------------------------------------
!
      MODULE tglf_coeff
!
! store the hermite basis matrix coefficients
!
      USE tglf_dimensions
      IMPLICIT NONE
      SAVE
! ave_h
      REAL :: hn,hp1,hp3,hr11,hr13,hr33
      REAl :: hw113,hw133,hw333
      REAL :: hu1,hu3,hu33,ht1,ht3
      REAL :: hu3ht1,hu3ht3,hu33ht1,hu33ht3
      REAL :: hb1,hb3,hb33,hb1ht1,hb3ht3,hb33ht1
      REAL :: hd1,hd3,hd33,hd1hu1,hd3hu3,hd33hu1
      REAL :: hv1r,hv2r,hv3r,hv4r,hv5r
      REAL :: hv6r,hv7r,hv8r,hv9r,hv10r
      REAL :: hv1rht1,hv2rht3,hv3rht1,hv4rht3
      REAL :: hv6rhu1,hv7rhu3,hv9rhu1,hv10rhu3
      REAL :: hv1i,hv2i,hv3i,hv4i,hv5i
      REAL :: hv6i,hv7i,hv8i,hv9i,hv10i
      REAL :: hv1iht1,hv2iht3,hv3iht1,hv4iht3
      REAL :: hv6ihu1,hv7ihu3,hv9ihu1,hv10ihu3
      REAL :: gradhp1,gradhr11,gradhr13
      REAL :: gradhp1p1,gradhr11p1,gradhr13p1
      REAl :: grad_hu1,grad_hu3
      REAL :: wdhu3ht1,wdhu3ht3,wdhu33ht1,wdhu33ht3
      REAL :: modwdhu3,modwdhu33
      REAL :: modwdhu3ht1,modwdhu3ht3,modwdhu33ht1,modwdhu33ht3
      REAL :: ave_hn(nsm,nb,nb)
      REAL :: ave_hp1(nsm,nb,nb),ave_hp3(nsm,nb,nb)
      REAL :: ave_hr11(nsm,nb,nb),ave_hr13(nsm,nb,nb),ave_hr33(nsm,nb,nb)
      REAL :: ave_hw113(nsm,nb,nb),ave_hw133(nsm,nb,nb),ave_hw333(nsm,nb,nb)
      REAL :: ave_ht1(nsm,nb,nb),ave_ht3(nsm,nb,nb)
      REAL :: ave_hu1(nsm,nb,nb),ave_hu3(nsm,nb,nb),ave_hu33(nsm,nb,nb)
      REAL :: ave_hu3ht1(nsm,nb,nb),ave_hu3ht3(nsm,nb,nb)
      REAL :: ave_hu33ht1(nsm,nb,nb),ave_hu33ht3(nsm,nb,nb)
      REAL :: ave_hninv(nsm,nb,nb),ave_hp1inv(nsm,nb,nb)
      REAL :: ave_hp3inv(nsm,nb,nb)
      REAL :: ave_gradhp1(nsm,nb,nb),ave_gradhr11(nsm,nb,nb)
      REAL :: ave_gradhr13(nsm,nb,nb)
      REAL :: ave_gradhp1p1(nsm,nb,nb),ave_gradhr11p1(nsm,nb,nb)
      REAL :: ave_gradhr13p1(nsm,nb,nb)
      REAL :: ave_gradhu1(nsm,nb,nb),ave_gradhu3(nsm,nb,nb)
      REAL :: ave_hnp0(nsm,nb,nb)
      REAL :: ave_hp1p0(nsm,nb,nb),ave_hp3p0(nsm,nb,nb)
      REAL :: ave_hr11p0(nsm,nb,nb),ave_hr13p0(nsm,nb,nb)
      REAL :: ave_hr33p0(nsm,nb,nb)
      REAL :: ave_hnb0(nsm,nb,nb)
      REAL :: ave_hp1b0(nsm,nb,nb),ave_hp3b0(nsm,nb,nb)
      REAL :: ave_hr11b0(nsm,nb,nb),ave_hr13b0(nsm,nb,nb)
      REAL :: ave_hr33b0(nsm,nb,nb)
      REAL :: ave_hw113b0(nsm,nb,nb),ave_hw133b0(nsm,nb,nb)
      REAl :: ave_hw333b0(nsm,nb,nb)
      REAL :: ave_gradhp1p0(nsm,nb,nb),ave_gradhr11p0(nsm,nb,nb)
      REAL :: ave_gradhr13p0(nsm,nb,nb)
      REAL :: ave_wdhu3ht1(nsm,nb,nb),ave_wdhu3ht3(nsm,nb,nb)
      REAL :: ave_wdhu33ht1(nsm,nb,nb),ave_wdhu33ht3(nsm,nb,nb)
      REAL :: ave_modwdhu3(nsm,nb,nb),ave_modwdhu33(nsm,nb,nb)
      REAL :: ave_modwdhu3ht1(nsm,nb,nb),ave_modwdhu3ht3(nsm,nb,nb)
      REAL :: ave_modwdhu33ht1(nsm,nb,nb),ave_modwdhu33ht3(nsm,nb,nb)
! ave_g
      REAL :: gn,gp1,gp3,gr11,gr13,gr33
      REAL :: gw113,gw133,gw333
      REAL :: gu1,gu3,gu33,gt1,gt3
      REAL :: gu3gt1,gu3gt3,gu33gt1,gu33gt3
      REAL :: gu1r,gu2r,gu3r,gu4r,gu5r
      REAL :: gu6r,gu7r,gu8r,gu9r,gu10r
      REAL :: gu1i,gu2i,gu3i,gu4i,gu5i
      REAL :: gu6i,gu7i,gu8i,gu9i,gu10i
      REAL :: gb1,gd1,gb3,gd3,gb33,gd33
      REAL :: gb1gt1,gd1gu1,gb3gt3,gd3gu3
      REAL :: gb33gt1,gd33gu1
      REAL :: gu1rgt1,gu1igt1,gu2rgt3,gu2igt3
      REAL :: gu3rgt1,gu3igt1,gu4rgt3,gu4igt3
      REAL :: gu6rgu1,gu7rgu3,gu9rgu1,gu10rgu3
      REAL :: gu6igu1,gu7igu3,gu9igu1,gu10igu3
      REAL :: gradgp1,gradgr11,gradgr13
      REAL :: grad_gu1,grad_gu3
      REAL :: gradgp1p1,gradgr11p1,gradgr13p1
      REAL :: wdgu3gt1,wdgu3gt3,wdgu33gt1,wdgu33gt3
      REAL :: modwdgu3,modwdgu33
      REAL :: modwdgu3gt1,modwdgu3gt3,modwdgu33gt1,modwdgu33gt3
      REAL :: ave_gn(nsm,nb,nb)
      REAL :: ave_gp1(nsm,nb,nb),ave_gp3(nsm,nb,nb)
      REAL :: ave_gr11(nsm,nb,nb),ave_gr13(nsm,nb,nb),ave_gr33(nsm,nb,nb)
      REAL :: ave_gw113(nsm,nb,nb),ave_gw133(nsm,nb,nb),ave_gw333(nsm,nb,nb)
      REAL :: ave_gt1(nsm,nb,nb),ave_gt3(nsm,nb,nb)
      REAL :: ave_gu1(nsm,nb,nb),ave_gu3(nsm,nb,nb),ave_gu33(nsm,nb,nb)
      REAL :: ave_gu3gt1(nsm,nb,nb),ave_gu3gt3(nsm,nb,nb)
      REAL :: ave_gu33gt1(nsm,nb,nb),ave_gu33gt3(nsm,nb,nb)
      REAL :: ave_gninv(nsm,nb,nb),ave_gp1inv(nsm,nb,nb)
      REAL :: ave_gp3inv(nsm,nb,nb)
      REAL :: ave_gradgp1(nsm,nb,nb),ave_gradgr11(nsm,nb,nb)
      REAL :: ave_gradgr13(nsm,nb,nb)
      REAL :: ave_gradgp1p1(nsm,nb,nb),ave_gradgr11p1(nsm,nb,nb)
      REAL :: ave_gradgr13p1(nsm,nb,nb)
      REAL :: ave_gradgu1(nsm,nb,nb),ave_gradgu3(nsm,nb,nb)
      REAL :: ave_gnp0(nsm,nb,nb)
      REAL :: ave_gp1p0(nsm,nb,nb),ave_gp3p0(nsm,nb,nb)
      REAL :: ave_gr11p0(nsm,nb,nb),ave_gr13p0(nsm,nb,nb)
      REAL :: ave_gr33p0(nsm,nb,nb)
      REAL :: ave_gnb0(nsm,nb,nb)
      REAL :: ave_gp1b0(nsm,nb,nb),ave_gp3b0(nsm,nb,nb)
      REAL :: ave_gr11b0(nsm,nb,nb),ave_gr13b0(nsm,nb,nb)
      REAL :: ave_gr33b0(nsm,nb,nb)
      REAL :: ave_gw113b0(nsm,nb,nb),ave_gw133b0(nsm,nb,nb)
      REAL :: ave_gw333b0(nsm,nb,nb)
      REAL :: ave_gradgp1p0(nsm,nb,nb),ave_gradgr11p0(nsm,nb,nb)
      REAL :: ave_gradgr13p0(nsm,nb,nb)
      REAL :: ave_wdgu3gt1(nsm,nb,nb),ave_wdgu3gt3(nsm,nb,nb)
      REAL :: ave_wdgu33gt1(nsm,nb,nb),ave_wdgu33gt3(nsm,nb,nb)
      REAL :: ave_modwdgu3(nsm,nb,nb),ave_modwdgu33(nsm,nb,nb)
      REAL :: ave_modwdgu3gt1(nsm,nb,nb),ave_modwdgu3gt3(nsm,nb,nb)
      REAL :: ave_modwdgu33gt1(nsm,nb,nb),ave_modwdgu33gt3(nsm,nb,nb)
! ave_wd_h
      REAL :: wdhp1,wdhp3,wdhr11,wdhr13,wdhr33
      REAL :: wdhu1,wdhu3,wdhu33,modwdhu1
      REAL :: wdht1,wdht3,modwdht1,modwdht3
      REAL :: ave_wdhp1(nsm,nb,nb),ave_wdhp3(nsm,nb,nb)
      REAL :: ave_wdhr11(nsm,nb,nb),ave_wdhr13(nsm,nb,nb)
      REAL :: ave_wdhr33(nsm,nb,nb)
      REAL :: ave_wdhu1(nsm,nb,nb),ave_wdhu3(nsm,nb,nb)
      REAL :: ave_wdhu33(nsm,nb,nb)
      REAL :: ave_modwdhu1(nsm,nb,nb)
      REAL :: ave_wdht1(nsm,nb,nb),ave_wdht3(nsm,nb,nb)
      REAL :: ave_modwdht1(nsm,nb,nb),ave_modwdht3(nsm,nb,nb)
! ave_wd_g
      REAL :: wdgp1,wdgp3,wdgr11,wdgr13,wdgr33
      REAL :: wdgu1,wdgu3,wdgu33,modwdgu1
      REAL :: wdgt1,wdgt3,modwdgt1,modwdgt3
      REAL :: ave_wdgp1(nsm,nb,nb),ave_wdgp3(nsm,nb,nb)
      REAL :: ave_wdgr11(nsm,nb,nb),ave_wdgr13(nsm,nb,nb)
      REAL :: ave_wdgr33(nsm,nb,nb)
      REAL :: ave_wdgu1(nsm,nb,nb),ave_wdgu3(nsm,nb,nb)
      REAL :: ave_wdgu33(nsm,nb,nb)
      REAL :: ave_modwdgu1(nsm,nb,nb)
      REAL :: ave_wdgt1(nsm,nb,nb),ave_wdgt3(nsm,nb,nb)
      REAL :: ave_modwdgt1(nsm,nb,nb),ave_modwdgt3(nsm,nb,nb)
! kpar_h
      COMPLEX :: kpar_hp1,kpar_hr11,kpar_hr13,kpar_hu1,kpar_hu3
      COMPLEX :: kpar_hb1,kpar_hb3,kpar_hb33
      COMPLEX :: kpar_hb1ht1,kpar_hb3ht3,kpar_hb33ht1
      REAL :: modkpar_hd1,modkpar_hd3,modkpar_hd33
      REAL :: modkpar_hd1hu1,modkpar_hd3hu3,modkpar_hd33hu1
      REAL :: modkpar_hu1,modkpar_hu3
      REAL :: ave_kparhp1(nsm,nb,nb),ave_kparhr11(nsm,nb,nb)
      REAL :: ave_kparhr13(nsm,nb,nb)
      REAL :: ave_kparhu1(nsm,nb,nb),ave_kparhu3(nsm,nb,nb)
      REAL :: ave_kparht1(nsm,nb,nb),ave_kparht3(nsm,nb,nb)
      REAL :: ave_modkparhu1(nsm,nb,nb),ave_modkparhu3(nsm,nb,nb)
! kpar_g
      COMPLEX :: kpar_gp1,kpar_gr11,kpar_gr13,kpar_gu1,kpar_gu3
      COMPLEX :: kpar_gb1,kpar_gb3,kpar_gb33
      COMPLEX :: kpar_gb1gt1,kpar_gb3gt3,kpar_gb33gt1
      REAL :: modkpar_gd1,modkpar_gd3,modkpar_gd33
      REAL :: modkpar_gd1gu1,modkpar_gd3gu3,modkpar_gd33gu1
      REAL :: modkpar_gu1,modkpar_gu3
      REAL :: ave_kpargp1(nsm,nb,nb),ave_kpargr11(nsm,nb,nb)
      REAL :: ave_kpargr13(nsm,nb,nb)
      REAL :: ave_kpargu1(nsm,nb,nb),ave_kpargu3(nsm,nb,nb)
      REAL :: ave_kpargt1(nsm,nb,nb),ave_kpargt3(nsm,nb,nb)
      REAL :: ave_modkpargu1(nsm,nb,nb),ave_modkpargu3(nsm,nb,nb)
! gradB_h
      REAL :: gradBhp1,gradBhp3
      REAL :: gradBhr11,gradBhr13,gradBhr33
      REAL :: gradBhu1,gradBhu3,gradBhu33
      REAL :: ave_gradBhp1(nsm,nb,nb),ave_gradBhp3(nsm,nb,nb)
      REAL :: ave_gradBhr11(nsm,nb,nb),ave_gradBhr13(nsm,nb,nb)
      REAL :: ave_gradBhr33(nsm,nb,nb)
      REAL :: ave_gradBhu1(nsm,nb,nb),ave_gradBhu3(nsm,nb,nb)
      REAL :: ave_gradBhu33(nsm,nb,nb)
! gradB_g
      REAL :: gradBgp1,gradBgp3
      REAL :: gradBgr11,gradBgr13,gradBgr33
      REAL :: gradBgu1,gradBgu3,gradBgu33
      REAL :: ave_gradBgp1(nsm,nb,nb),ave_gradBgp3(nsm,nb,nb)
      REAL :: ave_gradBgr11(nsm,nb,nb),ave_gradBgr13(nsm,nb,nb)
      REAL :: ave_gradBgr33(nsm,nb,nb)
      REAL :: ave_gradBgu1(nsm,nb,nb),ave_gradBgu3(nsm,nb,nb)
      REAL :: ave_gradBgu33(nsm,nb,nb)
!  ave_theta
      REAL :: gradB
      REAL :: ave_wd(nb,nb),ave_modwd(nb,nb)
      REAL :: ave_gradB(nb,nb),ave_lnB(nb,nb)
      REAL :: ave_b0(nb,nb),ave_b0inv(nb,nb)
      REAL :: ave_kpar(nb,nb),ave_modkpar(nb,nb)
      REAL :: ave_p0(nb,nb),ave_p0inv(nb,nb)
!
      END MODULE tglf_coeff
!
!-----------------------------------------------------
!
