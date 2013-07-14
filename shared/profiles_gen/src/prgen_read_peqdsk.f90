!--------------------------------------------------------------
! prgen_read_peqdsk.f90
!
! PURPOSE:
!  Extract peqdsk data from native TEXT file.
!--------------------------------------------------------------

subroutine prgen_read_peqdsk

  use prgen_globals

  implicit none

  character (len=100) :: t
  integer :: i
  real, dimension(:,:), allocatable :: xv
  integer, parameter :: ncol=3

  !----------------------------------------------------
  ! Read the peqdsk file
  !
  
  open(unit=1,file=raw_data_file,status='old')
  read(1,*) i, t
  peqdsk_nj = i
  nx = peqdsk_nj
  call allocate_peqdsk_vars

  ! psi_norm and ne(10^20/m^3)
  allocate(xv(ncol,i))
  read(1,*) xv
  peqdsk_psi(:) = xv(1,:)
  peqdsk_ne(:)  = xv(2,:)
  deallocate(xv)

  ! te(KeV)
  read(1,*) i, t
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_te,peqdsk_nj)
  deallocate(xv)

  ! ni(10^20/m^3)
  read(1,*) i, t
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_ni,peqdsk_nj)
  deallocate(xv)

  ! ti(KeV)
  read(1,*) i, t
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_ti,peqdsk_nj)
  deallocate(xv)

  if(peqdsk_ftype == 2) then

     ! nb(10^20/m^3)
     read(1,*) i, t
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_nb,peqdsk_nj)
     deallocate(xv)
     
     ! pb(KPa)
     read(1,*) i, t
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_pb,peqdsk_nj)
     deallocate(xv)
     
  endif

  ! omeg(kRad/s)
  read(1,*) i, t
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_omegat,peqdsk_nj)
  deallocate(xv)
  
  close(1)

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
