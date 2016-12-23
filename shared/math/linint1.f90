!------------------------------------------------------
! linint1.f90
!
! PURPOSE:
!  Given the gradient f of a function fi:
!
!                 x*
!                 /
!  f(x) = f(x*) + [ z(y) dy
!                 /
!                 x
!  
!  obtain f by trapezoidal integration.  It is assumed 
!  that f(x*) is known.
!
! NOTES:
!           z : (input) gradients
!           x : (input) radii
!          fb : (input) f(x*) = boundary condition 
!           f : (output) profile at x = f(x)
!------------------------------------------------------

subroutine linint1(z,x,n,fb,xs,f)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: z(n),x(n)
  real, intent(in) :: xs,fb
  real, intent(inout) :: f
  real :: x0
  
  integer :: i,i0
  real :: dx

  x0 = xs
  
  if (x0 > x(n)) then
     print *,'WARNING: (linint1) x0 > x(n)'
     x0 = x(n)
     f = fb
     return
  endif
  
  if (x0 < 0.0) then
     print *,'WARNING: (linint1) x0 < 0'
     x0 = 0.0
  endif
  
  ! Find i0 such that x(i0-1) < x0 < x(i0)
  do i=n,1,-1
     if (x0 > x(i)) exit
  enddo
  i0 = i+1
  if (i0 == 1) i0 = 2
  
  ! Integrate up to nearest gridpoint
  f = fb
  do i=n,i0+1,-1
     f = f-0.5*(z(i)+z(i-1))*(x(i)-x(i-1))
  enddo

  ! Final subinterval
  dx = x(i0)-x(i0-1)
  f = f-0.5/dx*(z(i0)*(dx**2-(x0-x(i0-1))**2)+z(i0-1)*(x0-x(i0))**2)
  
end subroutine linint1

