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
      integer  ms
      parameter(ms = 128)  ! ms needs to be even
! INPUT
      real*8 R(0:ms), Z(0:ms), Bp(0:ms)
      real*8 ds, Ls
      real*8 B_unit, Rmaj_s, rmin_s, q_s
      real*8 p_prime_s, q_prime_s
      real*8 p_prime_zero_s
! OUTPUT
      real*8 costheta_geo(0:ms),sintheta_geo(0:ms)
      real*8 costheta_p_geo(0:ms),s_p(0:ms)
      real*8 pk_geo(0:ms),epsl_geo(0:ms),qrat_geo(0:ms)
      real*8 kxoky_geo(0:ms), b_geo(0:ms),t_s(0:ms)
      real*8 S_prime(0:ms),kx_factor(0:ms)
      real*8 f,ff_prime
!
      common /s_grid/ R,Z,Bp,
     >  costheta_geo,sintheta_geo,costheta_p_geo,s_p,
     >  pk_geo,epsl_geo,qrat_geo,kxoky_geo,b_geo,t_s,
     >  S_prime,kx_factor,
     >  ds,Ls,B_unit,Rmaj_s, 
     >  rmin_s,q_s,p_prime_s,q_prime_s,p_prime_zero_s,
     >  f,ff_prime
!
! end of s_grid

  
