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

  !call prgen_get_chi(nx,q,dpsi,rho,torfluxa)
  
  !---------------------------------------------------------
  ! Map profile data into expro interface names:
  !
  expro_n_exp = nx
  expro_n_ion = 1
  call expro_init(1)

  expro_name(1) = 'D'
  expro_type(1) = '[therm]'
  expro_rho  = rho
  expro_rmin = rmin
  expro_rmaj = rmaj

 end subroutine prgen_map_null

