!----------------------------------------------------------------
! prgen_read_omfit.f90
!
! PURPOSE:
!  Read data from OMFIT/GACODE flux-surface mapper 
!----------------------------------------------------------------

subroutine prgen_read_omfit

  use prgen_globals

  implicit none

  integer :: i
  integer, parameter :: npsi=64
  real, dimension(3,npsi) :: efit_si,efit_ci,efit_xi,dummy3
  real, dimension(4,npsi) :: dummy4
  real, dimension(npsi) :: efit_rho,efit_psi,efit_q,efit_p
  real, dimension(npsi) :: efit_rmin,efit_rmaj,efit_kappa,efit_zmaj  
  real :: efit_torfluxa
  
  !----------------------------------------------------------------
  ! Read the fit files
  open(unit=1,file='bin.si.fit',status='old',access='stream')
  read(1) efit_si
  close(1)
  
  open(unit=1,file='bin.ci.fit',status='old',access='stream')
  read(1) efit_ci
  close(1)

  open(unit=1,file='bin.xi.fit',status='old',access='stream')
  read(1) dummy4
  close(1)
  
  open(unit=1,file='bin.psi.fit',status='old',access='stream')
  read(1) dummy3
  close(1)
  !----------------------------------------------------------------
  efit_psi   = dummy3(1,:)
  efit_q     = dummy3(2,:)
  efit_p     = dummy3(3,:)
  efit_rmin  = dummy4(1,:)
  efit_rmaj  = dummy4(2,:)
  efit_kappa = dummy4(3,:)/dummy4(1,:)
  efit_zmaj  = dummy4(4,:)
  
  call prgen_get_chi(npsi,efit_q,efit_psi,efit_rho,efit_torfluxa)
  call cub_spline(efit_rho,efit_psi,npsi,rho,dpsi,nx)
  call cub_spline(efit_rho,efit_rmin,npsi,rho,rmin,nx)

  call cub_spline(efit_psi,efit_q,npsi,dpsi,q,nx)
  call cub_spline(efit_psi,efit_rmaj,npsi,dpsi,rmaj,nx)
  call cub_spline(efit_psi,efit_kappa,npsi,dpsi,kappa,nx)
  call cub_spline(efit_psi,efit_si(1,:),npsi,dpsi,delta,nx) ; delta = sin(delta)
  call cub_spline(efit_psi,efit_si(2,:),npsi,dpsi,zeta,nx)

end subroutine prgen_read_omfit
