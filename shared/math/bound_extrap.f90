!-------------------------------------------------------
! bound_extrap.f90
!
! PURPOSE:
!  EXtrapolate to left and right boundaries. 
!-------------------------------------------------------

subroutine bound_extrap(fa,fb,f,r,n)

  implicit none

  integer, intent(in) :: n

  real, intent(inout) :: fa
  real, intent(inout) :: fb
  real, intent(in), dimension(n) :: f
  real, intent(in), dimension(n) :: r

  real :: r1,r2,r3
  real :: ra,rb
  real :: f1,f2,f3


  ! Left boundary

  ra = r(1)

  r2 = r(2)
  r3 = r(3)
  
  f2 = f(2)
  f3 = f(3)

  fa = (ra-r2)/(r3-r2)*f3+(r3-ra)/(r3-r2)*f2

  ! Right boundary

  rb = r(n)

  r1 = r(n-2)
  r2 = r(n-1)
  
  f1 = f(n-2)
  f2 = f(n-1)

  fb = (rb-r1)/(r2-r1)*f2+(r2-rb)/(r2-r1)*f1

end subroutine bound_extrap
