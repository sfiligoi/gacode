!------------------------------------------------------------------------------
! cgyro_fourier_extend.f90
!
! PURPOSE:
!  Compute coefficients of complex DFT over interval x=[0,2pi]
! 
!            inf        ilx
!    f(x) =  Sum  c(l) e   
!           l=-inf           
!
! Given f(x), compute c(l).
!
! EXAMPLE: 
! 
!             inf   2                      ilx              1
!     x-pi =  Sum  --- sin(lx) = Sum c(l) e    with c(l) = ---
!             l=1   l             l                         il
!------------------------------------------------------------------------------

subroutine cgyro_fourier_extend

  use cgyro_globals

  implicit none

  integer,parameter :: nxi=512
  integer :: l,i
  real, dimension(:), allocatable :: x,f
  real, external :: gl_f,gl_fp
  real :: u,a,b,c,d
  real :: e
  real :: f0,g0,f1,g1

  allocate(x(nxi),f(nxi))

  do i=1,nxi
     x(i) = 2*pi*(i-1.0)/nxi
  enddo

  ! Scale fraction to domain (0,2pi)
  e = eps_global*2*pi

  ! Function and derivative conditions
  f0 = gl_f(e)
  f1 = gl_fp(e)
  g0 = gl_f(2*pi-e)
  g1 = gl_fp(2*pi-e)

  ! Polynomial coefficients
  c = (f1-g1)/(4*e)
  a = (f0+g0)/2-c*e**2
  d = (e*(f1+g1)-(f0-g0))/(4*e**3)
  b = ((f0-g0)-2*d*e**3)/(2*e)

  ! Definition of extended function 
  do i=1,nxi
     if (x(i) <= e) then
        u = x(i)
        f(i) = a+b*u+c*u**2+d*u**3
     else if (x(i) >= 2*pi-e) then
        u = x(i)-2*pi
        f(i) = a+b*u+c*u**2+d*u**3
     else
        f(i) = gl_f(x(i))
     endif
  enddo

  ! Project complex Fourier coefficients cc(l)
  do l=-n_global,n_global
     cg(l) = sum(f(:)*exp(i_c*l*x(:)))/nxi
  enddo
  deallocate(x,f)

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_extend,status='replace')
     write(io,*) 'Extended-domain Fourier series'
     write(io,*) 
     write(io,*) 'ExB shear: c_l ~ 1/(il)'
     write(io,*)
     write(io,'(a)') '  l     Re(c_l)       Im(c_l)' 
     do l=-n_global,n_global
        write(io,'(i3,2x,2(1pe13.6,1x))') l,cg(l)
     enddo
     close(io)
  endif

end subroutine cgyro_fourier_extend

real function gl_f(x)

  real, intent(in) :: x
  real, parameter :: pi = 3.1415926535897932

  gl_f = x-pi

end function gl_f

real function gl_fp(x)

  real, intent(in) :: x

  gl_fp = 1.0

end function gl_fp


