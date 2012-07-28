!--------------------------------------------------------------
! prgen_read_ufile.f90
!
! PURPOSE:
!  Extract data from UFILE using external library (read
!--------------------------------------------------------------

subroutine prgen_read_ufile

  use prgen_read_globals

  implicit none

  integer :: i
  real :: chita

  open(unit=1,file='out.BT.ave',status='old') 
  read(1,*) ufile_bref
  close(1)
  open(unit=1,file='out.PHIA.ave',status='old') 
  read(1,*) chita
  close(1)
  open(unit=1,file='out.com',status='old') 
  read(1,*) raw_data_file
  read(1,*) ufile_tok
  read(1,*) ufile_shot
  close(1)

  ! Use classic parameterization chi_t = B_ref/2 rho^2
  !
  ufile_arho = sqrt(2.0*chita/ufile_bref)

  nx = 6

  call allocate_internals
  call allocate_ufile_vars

  do i=1,nx
     rho(i) = (i-1)/(nx-1.0)
  enddo
 
  call ufile_mapper('out.NE.ave',rho,ufile_ne,nx)
  call ufile_mapper('out.TE.ave',rho,ufile_te,nx)
  call ufile_mapper('out.ZEFFR.ave',rho,ufile_zeff,nx)
  call ufile_mapper('out.RMAJOR.ave',rho,rmaj,nx)
  call ufile_mapper('out.RMINOR.ave',rho,rmin,nx)
  call ufile_mapper('out.Q.ave',rho,q,nx)
  call ufile_mapper('out.KAPPAR.ave',rho,kappa,nx)
  call ufile_mapper('out.DELTAR.ave',rho,delta,nx)

end subroutine prgen_read_ufile

subroutine ufile_mapper(file,x,y,nx)

  character(len=*), intent(in) :: file
  integer, intent(in) :: nx
  real, intent(in) :: x(nx)
  real, intent(inout) :: y(nx)

  real, dimension(:), allocatable :: x0
  real, dimension(:), allocatable :: y0
  real :: ya,yb
  integer :: nx0

  open(unit=1,file='out.dim',status='old')
  read(1,*) nx0
  close(1)

  allocate(x0(0:nx0+1))
  allocate(y0(0:nx0+1))

  open(unit=1,file=trim(file),status='old')
  do i=1,nx0
     read(1,*) x0(i),y0(i)
  enddo

  if (x0(1) > 0.0 .and. x0(nx0) < 1.0) then
     ! Extrapolate to 0 and 1
     x0(0)     = 0.0
     x0(nx0+1) = 1.0
     call bound_extrap(ya,yb,y0,x0,nx0+2)
     y0(0) = ya
     y0(nx0+1) = yb
     call cub_spline(x0,y0,nx0+2,x,y,nx)
  else if (x0(1) > 0.0 .and. x0(nx0) == 1.0) then
     ! Extrapolate to 0 only
     x0(0) = 0.0
     call bound_extrap(ya,yb,y0(0:nx0),x0(0:nx0),nx0+1)
     y0(0) = ya
     call cub_spline(x0(0:nx0),y0(0:nx0),nx0+1,x,y,nx)
  else
     print *,'ERROR IN UFILE_MAPPER'
     stop
  endif

  deallocate(x0)
  deallocate(y0)

end subroutine ufile_mapper
