!--------------------------------------------------------------
! prgen_read_peqdsk.f90
!
! PURPOSE:
!  Extract peqdsk data from native TEXT file.
!--------------------------------------------------------------

subroutine prgen_read_peqdsk

  use prgen_globals

  implicit none

  integer, parameter :: ncol=3
  integer :: i
  logical :: ierr
  real, dimension(:,:), allocatable :: xv

  !----------------------------------------------------
  ! Read the peqdsk file
  !

  open(unit=1,file='pfile.ne',status='old')
  read(1,*) i
  peqdsk_nj = i
  nx = peqdsk_nj
  call allocate_peqdsk_vars
  close(1)

  ! psi_norm and ne(10^20/m^3)
  open(unit=1,file='pfile.ne',status='old')
  allocate(xv(ncol,i))
  read(1,*) i
  read(1,*) xv
  peqdsk_psi(:) = xv(1,:)
  peqdsk_ne(:)  = xv(2,:)
  deallocate(xv)
  close(1)

  ! te(KeV)
  open(unit=1,file='pfile.te',status='old')
  read(1,*) i
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_te,peqdsk_nj)
  deallocate(xv)
  close(1)

  ! ni(10^20/m^3)
  open(unit=1,file='pfile.ni',status='old')
  read(1,*) i
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_ni,peqdsk_nj)
  deallocate(xv)
  close(1)

  ! ti(KeV)
  open(unit=1,file='pfile.ti',status='old')
  read(1,*) i
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_ti,peqdsk_nj)
  deallocate(xv)
  close(1)

  ! nb(10^20/m^3)
  inquire(file='pfile.nb',exist=ierr)
  if(ierr) then
     open(unit=1,file='pfile.nb',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_nb,peqdsk_nj)
     deallocate(xv)
     close(1)
  else
     peqdsk_nb(:) = 0.0
  endif

  ! pb(KPa)
  inquire(file='pfile.pb',exist=ierr)
  if(ierr) then
     open(unit=1,file='pfile.pb',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_pb,peqdsk_nj)
     deallocate(xv)
     close(1)
  else
     peqdsk_pb(:) = 0.0
  endif

  ! omeg(kRad/s)
  inquire(file='pfile.omeg',exist=ierr)
  if(ierr) then
     open(unit=1,file='pfile.omeg',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_omegat,peqdsk_nj)
     deallocate(xv)
     close(1)
  else
     peqdsk_omegat(:) = 0.0
  endif

  ! omgeb(kRad/s)
  inquire(file='pfile.omgeb',exist=ierr)
  if(ierr) then
     open(unit=1,file='pfile.omgeb',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_omgeb,peqdsk_nj)
     deallocate(xv)
     close(1)
  else
     peqdsk_omgeb(:) = 0.0
  endif

  call allocate_internals

  dpsi(:)   = peqdsk_psi(:)-peqdsk_psi(1)
  rmin(:)   = 0.0
  rmaj(:)   = 0.0
  rho(:)    = 0.0
  kappa(:)  = 0.0
  delta(:)  = 0.0
  zmag(:)   = 0.0
  zeta(:)   = 0.0
  omega0(:) = 0.0

end subroutine prgen_read_peqdsk
