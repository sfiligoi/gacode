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
  integer :: npsi
  real, dimension(:,:), allocatable :: efit_si,efit_ci
  real, dimension(:), allocatable :: efit_rho,efit_psi,efit_q,efit_p
  real, dimension(:), allocatable :: efit_rmin,efit_rmaj,efit_kappa,efit_zmaj  
  real :: efit_torfluxa
  
  !----------------------------------------------------------------
  open(unit=1,file='out.dim',status='old')
  read(1,*) npsi
  close(1)

  allocate(efit_si(3,npsi))
  allocate(efit_ci(3,npsi))
  allocate(efit_rmin(npsi))
  allocate(efit_rmaj(npsi))
  allocate(efit_kappa(npsi))
  allocate(efit_zmaj(npsi))
  allocate(efit_psi(npsi))
  allocate(efit_q(npsi))
  allocate(efit_p(npsi))
  allocate(efit_rho(npsi))

  open(unit=1,file='out.data',status='old',access='stream')
  read(1) efit_si
  read(1) efit_ci
  read(1) efit_rmin
  read(1) efit_rmaj
  read(1) efit_kappa
  read(1) efit_zmaj
  read(1) efit_psi
  read(1) efit_q
  read(1) efit_p
  close(1)
  !----------------------------------------------------------------
  
  efit_psi  = efit_psi-efit_psi(1)
  dpsi_data = dpsi(nx)
  dpsi_efit = efit_psi(npsi)
  
  call prgen_get_chi(npsi,efit_q,efit_psi,efit_rho,efit_torfluxa)
  call cub_spline(efit_rho,efit_psi,npsi,rho,dpsi,nx)
  call cub_spline(efit_rho,efit_rmin,npsi,rho,rmin,nx)

  call cub_spline(efit_psi,efit_q,npsi,dpsi,q,nx)
  call cub_spline(efit_psi,efit_p,npsi,dpsi,p_tot,nx)
  call cub_spline(efit_psi,efit_zmaj,npsi,dpsi,zmag,nx)
  call cub_spline(efit_psi,efit_rmaj,npsi,dpsi,rmaj,nx)
  call cub_spline(efit_psi,efit_kappa,npsi,dpsi,kappa,nx)
  call cub_spline(efit_psi,efit_si(2,:),npsi,dpsi,delta,nx) ; delta = sin(delta)
  call cub_spline(efit_psi,efit_si(3,:),npsi,dpsi,zeta,nx)

end subroutine prgen_read_omfit
