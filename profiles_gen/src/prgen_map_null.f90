!------------------------------------------------------------
! prgen_map_null.f90
!
! PURPOSE:
!  Map shape data only.  
!------------------------------------------------------------

subroutine prgen_map_null

  use prgen_globals
  use expro

  implicit none

  call prgen_get_chi(nx,q,dpsi,rho,expro_torfluxa)

  !---------------------------------------------------------
  ! Map profile data into expro interface names:
  !
  expro_n_exp = nx
  expro_n_ion = 1
  call expro_init(1)
  
  expro_rho  = rho
  expro_rmin = rmin(:)
  expro_rmaj = rmaj(:)
  ! COORDINATES: set sign of q
  expro_q = abs(q(:))*ipccw*btccw
  expro_ptot = p_tot(:)
  expro_kappa = kappa(:)
  expro_delta = delta(:)
  expro_zeta = zeta(:)
  expro_zmag = zmag(:)
  expro_te = 1.0
  expro_ne = 1.0

  ! COORDINATES: set sign of poloidal flux
  expro_polflux = abs(dpsi(:))*(-ipccw)

  ! Ion temperatures and densities
  expro_ni(1,:) = 1.0
  expro_ti(1,:) = 1.0

 end subroutine prgen_map_null

