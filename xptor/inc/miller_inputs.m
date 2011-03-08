!---------------------------------------------------------------
! MILLER_inputs.m
!
! PURPOSE:
!  include file which contains input parameters passed 
!  though a call to miller_init.f
!
! REVISIONS:
! 24 June 05: gms
!  this version for TGLF separated miller and mercier-luc components
!  from GKS and GYRO versions
!---------------------------------------------------------------
      real*8 rmin_loc
      real*8 rmaj_loc
      real*8 aspectratio_loc
      real*8 q_loc
      real*8 delta_loc
      real*8 kappa_loc
      real*8 shift_loc
      real*8 shat_loc
      real*8 beta_loc
      real*8 s_delta_loc
      real*8 s_kappa_loc
      real*8 dlntidr_loc
      real*8 dlntedr_loc
      real*8 dlnnidr_loc
      real*8 dlnnedr_loc
      real*8 tiote_loc
      real*8 nione_loc
      real*8 shat_mhd_loc
      real*8 alpha_mhd_loc
      real*8 volume_loc
      real*8 volume_prime
! Flags
      real*8 dlnpdr_loc
      real*8 p_prime_zero
      integer alpha_mhd_method
      integer curve_method
      integer alpha_flag
      integer pressure_flag
!
      common /miller_in/ rmin_loc, rmaj_loc, aspectratio_loc, q_loc,
     > delta_loc, kappa_loc, shift_loc, shat_loc, beta_loc, s_delta_loc,
     > s_kappa_loc, dlntidr_loc, dlntedr_loc, dlnnidr_loc, dlnnedr_loc, 
     > tiote_loc, nione_loc, shat_mhd_loc, alpha_mhd_loc, dlnpdr_loc,
     > p_prime_zero, alpha_mhd_method, curve_method,
     > alpha_flag, pressure_flag 
!
  
