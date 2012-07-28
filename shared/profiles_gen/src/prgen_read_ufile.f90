!--------------------------------------------------------------
! prgen_read_ufile.f90
!
! PURPOSE:
!  Read profiles from post-processed (by gacode_ufile_tool) 
!  UFILEs. 
!--------------------------------------------------------------

subroutine prgen_read_ufile

  use prgen_read_globals

  implicit none

  integer :: i
  real :: chi_ta
  real, dimension(:), allocatable :: chi_t

  open(unit=1,file='out.BT.ave',status='old') 
  read(1,*) ufile_bref
  close(1)

  ! Chi_t(a), where Chi_t is the total toroidal flux
  open(unit=1,file='out.PHIA.ave',status='old') 
  read(1,*) chi_ta
  chi_ta = chi_ta/(2*pi)
  close(1)

  open(unit=1,file='out.com',status='old') 
  read(1,*) raw_data_file
  read(1,*) ufile_tok
  read(1,*) ufile_shot
  close(1)

  nx = n_ufile

  call allocate_internals
  call allocate_ufile_vars

  allocate(chi_t(nx))

  ! Classic normalized rho
  do i=1,nx
     rho(i) = (i-1)/(nx-1.0)
  enddo
 
  call ufile_mapper('out.NE.ave',rho,ufile_ne,nx,1)
  call ufile_mapper('out.TE.ave',rho,ufile_te,nx,1)
  call ufile_mapper('out.TI.ave',rho,ufile_ti,nx,1)
  call ufile_mapper('out.ZEFFR.ave',rho,ufile_zeff,nx,0)
  call ufile_mapper('out.RMAJOR.ave',rho,rmaj,nx,0)
  call ufile_mapper('out.RMINOR.ave',rho,rmin,nx,0)
  call ufile_mapper('out.Q.ave',rho,q,nx,0)
  call ufile_mapper('out.KAPPAR.ave',rho,kappa,nx,0)
  call ufile_mapper('out.DELTAR.ave',rho,delta,nx,0)
  call ufile_mapper('out.PRES.ave',rho,ufile_pres,nx,1)

  call ufile_mapper('out.NM1.ave',rho,ufile_nm1,nx,1)
  ufile_nion=1
  call ufile_mapper('out.NM2.ave',rho,ufile_nm2,nx,1)
  if (ufile_nm2(1) > 0.0) ufile_nion=2 
  call ufile_mapper('out.NM3.ave',rho,ufile_nm3,nx,1)
  if (ufile_nm3(1) > 0.0) ufile_nion=3
  
  !------------------------------------------------------------------------
  ! Use classic parameterization chi_t = B_ref/2 rho^2
  ! where chi_t is toroidal flux over 2pi.
  !
  ufile_arho = sqrt(2*chi_ta/ufile_bref)
  chi_t(:)   = chi_ta*rho**2
  !
  ! Compute psi_p by integrating d(chi_t)/q = d(psi_p)
  ! using the trapezoidal rule
  !
  dpsi(1) = 0.0
  do i=2,nx
     dpsi(i) = dpsi(i-1) + 0.5*(1.0/q(i-1)+1.0/q(i))*(chi_t(i)-chi_t(i-1))
  enddo
  !------------------------------------------------------------------------

end subroutine prgen_read_ufile

subroutine ufile_mapper(file,x,y,nx,neg_protect)

  character(len=*), intent(in) :: file
  integer, intent(in) :: nx
  integer, intent(in) :: neg_protect
  real, intent(in) :: x(nx)
  real, intent(inout) :: y(nx)

  real, dimension(:), allocatable :: x0
  real, dimension(:), allocatable :: y0
  real :: ya,yb
  integer :: nx0
  integer :: ierr

  open(unit=1,file='out.dim',status='old')
  read(1,*) nx0
  close(1)

  allocate(x0(0:nx0+1))
  allocate(y0(0:nx0+1))

  open(unit=1,file=trim(file),status='old',iostat=ierr)
  if (ierr == 0) then
     ! File exists
     do i=1,nx0
        read(1,*) x0(i),y0(i)
     enddo
  else
     ! File does NOT exist: pad output with zeros
     y(:) = 0.0
     return
  endif

  if (x0(1) > 0.0 .and. x0(nx0) < 1.0) then
     ! Extrapolate to 0 and 1
     x0(0)     = 0.0
     x0(nx0+1) = 1.0
     if (neg_protect==0) then
        call bound_extrap(ya,yb,y0,x0,nx0+2)
     else
        call bound_exp(ya,yb,y0,x0,nx0+2)
     endif
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
     print *,'ERROR: (prgen) Boundary issue in prgen_read_ufile.'
     stop
  endif

  deallocate(x0)
  deallocate(y0)

end subroutine ufile_mapper

subroutine bound_exp(fa,fb,f,r,n)

  implicit none

  integer, intent(in) :: n

  real, intent(inout) :: fa
  real, intent(inout) :: fb
  real, intent(in), dimension(n) :: f
  real, intent(in), dimension(n) :: r

  real :: r1,r2,r3
  real :: ra,rb
  real :: f1,f2,f3
  real :: l,c

  ! Left boundary

  ra = r(1)

  r2 = r(2)
  r3 = r(3)
  
  f2 = f(2)
  f3 = f(3)

  fa = (ra-r2)/(r3-r2)*f3+(r3-ra)/(r3-r2)*f2

  ! Right boundary (fit to exp(-r/l) with scale length l)

  rb = r(n)

  r1 = r(n-2)
  r2 = r(n-1)
  
  f1 = f(n-2)
  f2 = f(n-1)

  l = (r2-r1)/log(f1/f2)
  c = f2*exp(r2/l)

  fb = c*exp(-rb/l)

end subroutine bound_exp
