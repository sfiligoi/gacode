!------------------------------------------------
! EXPRO_interface.f90
!
! PURPOSE:
!  Interface to experimental profiles.
!
! NOTES:
!  For more detailed variable definitions, see
!
!  http://fusion.gat.com/theory/input.profiles
!  
!  * Scalars:
!
!  EXPRO_ncol
!  EXPRO_nblock
!  EXPRO_n_exp
!  EXPRO_b_ref
!  EXPRO_arho
!
!  * Vectors:
!
!  EXPRO_rho(:)
!  EXPRO_rmin(:)
!  EXPRO_rmaj(:)
!  EXPRO_q(:)
!  EXPRO_kappa(:)
!
!  EXPRO_delta(:)
!  EXPRO_te(:)
!  EXPRO_ne(:)
!  EXPRO_z_eff(:)
!  EXPRO_w0(:)
!
!  EXPRO_flow_mom(:)
!  EXPRO_pow_e(:)
!  EXPRO_pow_i(:)
!  EXPRO_pow_ei(:)
!  EXPRO_zeta(:)
!
!  EXPRO_flow_beam(:)
!  EXPRO_flow_wall(:)
!  EXPRO_zmag(:)
!  EXPRO_ptot(:) ! WG
!  [EMPTY]
!
!  EXPRO_ni(1,:)
!  EXPRO_ni(2,:)
!  EXPRO_ni(3,:)
!  EXPRO_ni(4,:)
!  EXPRO_ni(5,:)
!
!  EXPRO_ti(1,:)
!  EXPRO_ti(2,:)
!  EXPRO_ti(3,:)
!  EXPRO_ti(4,:)
!  EXPRO_ti(5,:)
!
!  EXPRO_vtor(1,:)
!  EXPRO_vtor(2,:)
!  EXPRO_vtor(3,:)
!  EXPRO_vtor(4,:)
!  EXPRO_vtor(5,:)
!
!  EXPRO_vpol(1,:)
!  EXPRO_vpol(2,:)
!  EXPRO_vpol(3,:)
!  EXPRO_vpol(4,:)
!  EXPRO_vpol(5,:)
!
!  * Control parameters (user can change these)
! 
!  EXPRO_ctrl_density_method (1=do nothing, 2=force quasin.)
!  EXPRO_ctrl_z(1:5) (ion charges)
!  EXPRO_ctrl_numeq_flag (0=model,1=numerical)
!  EXPRO_ctrl_signq 
!  EXPRO_ctrl_signb
!  EXPRO_ctrl_rotation_method (1=candy-phi,2=waltz-U_parallel)
!
!  * Derived quantities:
!
!  EXPRO_bunit(:)       B_unit (T)
!  EXPRO_s(:)           (r/q)(dq/dr)
!  EXPRO_drmaj(:)       dR_0/dr
!  EXPRO_dzmag(:)       dZ_0/dr
!  EXPRO_sdelta(:)      r d(delta)/dr
!  EXPRO_skappa(:)      (r/kappa) d(kappa)/dr
!  EXPRO_szeta(:)       r d(zeta)/dr
!  EXPRO_dlnnedr(:)     -dln(ne)/dr (1/m)
!  EXPRO_dlntedr(:)     -dln(Te)/dr (1/m)
!  EXPRO_dlnnidr(1:5,:) -dln(ni)/dr (1/m)
!  EXPRO_dlntidr(1:5,:) -dln(ti)/dr (1/m)
!  EXPR0_dlnptotdr(:)   -dln(ptot)/dr (1/m) ! WG
!  EXPRO_w0(:)          w0 (1/s)
!  EXPRO_w0p(:)         d(w0)/dr (1/s/m)
!  EXPRO_vol(:)         V (m^3)
!  EXPRO_volp(:)        dV/dr (m^2)
!  EXPRO_cs(:)          cs (m/s)
!  EXPRO_rhos(:)        rhos (m)
!  EXPRO_ni_new(:)      ni [Corrected for quasin.]
!  EXPRO_dlnnidr_new(:) -dln(ni)/dr (1/m) [Corrected for quasin.] 
!  EXPRO_grad_r0(:)     |grad r| at theta=0 
!  EXPRO_bp0(:)         B_pol at theta=0 (T)
!  EXPRO_bt0(:)         B_tor at theta=0 (T)
!  EXPRO_poloidalfluxover2pi(:) psi_p (T*m^2/rad=Wb/rad)
!  EXPRO_gamma_e(:)     r/q d(w0)/dr (1/s)
!  EXPRO_gamma_p(:)     R_0 d(w0)/dr (1/s)
!  EXPRO_mach(:)        R_0 w0/cs
!
!  * General geometry arrays:
!
!  EXPRO_geo(:,:,:)
!  EXPRO_dgeo(:,:,:)
!------------------------------------------------ 

module EXPRO_interface

  integer :: EXPRO_error=0
  integer, parameter :: nion_max=5

  ! Fundamental input.profiles scalars

  integer :: EXPRO_ncol
  integer :: EXPRO_nblock
  integer :: EXPRO_n_exp
  real    :: EXPRO_b_ref
  real    :: EXPRO_arho

  ! Fundamental input.profiles arrays

  real, dimension(:),allocatable :: EXPRO_rho
  real, dimension(:),allocatable :: EXPRO_rmin
  real, dimension(:),allocatable :: EXPRO_rmaj
  real, dimension(:),allocatable :: EXPRO_q
  real, dimension(:),allocatable :: EXPRO_kappa

  real, dimension(:),allocatable :: EXPRO_delta
  real, dimension(:),allocatable :: EXPRO_te
  real, dimension(:),allocatable :: EXPRO_ne
  real, dimension(:),allocatable :: EXPRO_z_eff
  real, dimension(:),allocatable :: EXPRO_w0

  real, dimension(:),allocatable :: EXPRO_flow_mom
  real, dimension(:),allocatable :: EXPRO_pow_e
  real, dimension(:),allocatable :: EXPRO_pow_i
  real, dimension(:),allocatable :: EXPRO_pow_ei
  real, dimension(:),allocatable :: EXPRO_zeta
 
  real, dimension(:),allocatable :: EXPRO_flow_beam
  real, dimension(:),allocatable :: EXPRO_flow_wall
  real, dimension(:),allocatable :: EXPRO_zmag
  real, dimension(:),allocatable :: EXPRO_ptot ! WG

  real, dimension(:,:),allocatable :: EXPRO_ni
  real, dimension(:,:),allocatable :: EXPRO_ti
  real, dimension(:,:),allocatable :: EXPRO_vtor
  real, dimension(:,:),allocatable :: EXPRO_vpol

  ! Derived quantities

  real, dimension(:),allocatable :: EXPRO_bunit
  real, dimension(:),allocatable :: EXPRO_s
  real, dimension(:),allocatable :: EXPRO_drmaj
  real, dimension(:),allocatable :: EXPRO_dzmag
  real, dimension(:),allocatable :: EXPRO_sdelta
  real, dimension(:),allocatable :: EXPRO_skappa
  real, dimension(:),allocatable :: EXPRO_szeta
  real, dimension(:),allocatable :: EXPRO_dlnnedr
  real, dimension(:),allocatable :: EXPRO_dlntedr
  real, dimension(:,:),allocatable :: EXPRO_dlnnidr
  real, dimension(:,:),allocatable :: EXPRO_dlntidr
  real, dimension(:),allocatable :: EXPRO_dlnptotdr ! WG

  real, dimension(:),allocatable :: EXPRO_w0p

  real, dimension(:),allocatable :: EXPRO_vol
  real, dimension(:),allocatable :: EXPRO_volp

  real, dimension(:),allocatable :: EXPRO_cs
  real, dimension(:),allocatable :: EXPRO_rhos

  real, dimension(:),allocatable :: EXPRO_ni_new
  real, dimension(:),allocatable :: EXPRO_dlnnidr_new
  real, dimension(:),allocatable :: EXPRO_grad_r0

  real, dimension(:),allocatable :: EXPRO_bp0
  real, dimension(:),allocatable :: EXPRO_bt0
  real, dimension(:),allocatable :: EXPRO_poloidalfluxover2pi

  ! input.profiles.geo dimension and arrays

  integer :: EXPRO_nfourier
  real, dimension(:,:,:),allocatable :: EXPRO_geo
  real, dimension(:,:,:),allocatable :: EXPRO_dgeo

  ! Gyrokinetic rotation parameters

  real, dimension(:),allocatable :: EXPRO_gamma_e
  real, dimension(:),allocatable :: EXPRO_gamma_p
  real, dimension(:),allocatable :: EXPRO_mach
  
  ! Control parameters (force nonsensical defaults for usage check)
  
  integer :: EXPRO_ctrl_density_method = -1
  real, dimension(5) :: EXPRO_ctrl_z
  integer :: EXPRO_ctrl_numeq_flag=-1
  real :: EXPRO_ctrl_signb = -10.0
  real :: EXPRO_ctrl_signq = -10.0
  integer :: EXPRO_ctrl_rotation_method = -1
  
end module EXPRO_interface
