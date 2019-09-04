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
  integer :: npsi=64
  real, dimension(3,64) :: data_si,data_ci
  real, dimension(3,64) :: data_psi
  
  !----------------------------------------------------------------
  ! Read the OMFIT fit file
  open(unit=1,file='bin.si.fit',status='old',access='stream')
  read(1) data_si
  close(1)
  open(unit=1,file='bin.ci.fit',status='old',access='stream')
  read(1) data_ci
  close(1)
  open(unit=1,file='bin.psi.fit',status='old',access='stream')
  read(1) data_psi
  close(1)

  # delta
  call cub_spline(data_psi(0,:),data_si(1,:),npsi,dpsi,delta,nx)

  


end subroutine prgen_read_omfit
