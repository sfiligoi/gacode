!--------------------------------------------------------------------
! tgyro_scaleint.f90
!
! PURPOSE:
!  Given the linear (z=f') or negative logarithmic (z=-f'/f) gradient 
!  of a function f, find the backward integral:
!
! type=log
!                    x*
!                    /
!  f(x) = f(x*) exp[ | z(s) ds ]
!                    /
!                    x
!
! type=lin
!                 x*
!                 /
!  f(x) = f(x*) - | z(s) ds
!                 /
!                 x
!
!
!  In other words, We obtain f by analytic integration of the piecewise
!  linear gradient function z(x).
!
! NOTES:
!        z(n) : (input) gradients
!        x(n) : (input) radii
!           n : (input) dimension of input arrays
!          fb : (input) f(x*) = boundary condition 
!          x0 : (input) desired location of profile
!           f : (output) f(x0)
!---------------------------------------------------------------------

subroutine math_scaleint(z,x,n,fb,x0,f,type)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: z(n),x(n)
  real, intent(in) :: x0,fb
  real, intent(inout) :: f
  character(len=3), intent(in) :: type

  integer :: i,i0
  real :: xs,dx

  xs = x0

  if (xs > x(n)) then
     print *,'WARNING: (tgyro_scaleint) x0 > x(n)'
     xs = x(n)
     f = fb
     return
  endif

  if (xs < 0.0) then
     print *,'WARNING: (tgyro_scaleint) x0 < 0'
     xs = 0.0
  endif

  ! Find i0 such that x(i0-1) < xs < x(i0)
  do i0=n,2,-1
     if (xs >= x(i0-1)) exit
  enddo

  ! Integrate up to nearest gridpoint with trapezoidal integration
  ! then perform final subinterval with analytic integration.
  ! NOTE: the analytic method reduces to trapezoidal rule at meshpoints

  dx = x(i0)-x(i0-1)
  f  = fb

  if (type == 'log') then
     do i=n,i0+1,-1
        f = f*exp(0.5*(z(i)+z(i-1))*(x(i)-x(i-1)))
     enddo
     ! Analytic integral of linear function
     f = f*exp(0.5/dx*(z(i0)*(dx**2-(xs-x(i0-1))**2)+z(i0-1)*(xs-x(i0))**2))
  else
     do i=n,i0+1,-1
        f = f-0.5*(z(i)+z(i-1))*(x(i)-x(i-1))
     enddo
     ! Analytic integral of linear function
     f = f-0.5/dx*(z(i0)*(dx**2-(xs-x(i0-1))**2)+z(i0-1)*(xs-x(i0))**2)
  endif

end subroutine math_scaleint

