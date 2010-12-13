!------------------------------------------------
! invert_p32.f90
!
! PURPOSE:
!  Calculate x given s, where:
!
!                         x
!                  2      /          -t
! s = p32(x) = --------   | sqrt(t) e   dt
!              sqrt(pi)   /
!                         0
!
! Speed is not esential, so we use bisection 
! for robsutness.
!
! REVISIONS:
! 26 Sept 01: jc 
!  Documented
! 25 May 05: jc
!  Moved into gyro/math subdir and renamed.
!------------------------------------------------

subroutine invert_p32(s,x,x_max)

  implicit none

  real, parameter :: eps_x = 1e-10

  real, intent(in) :: s
  real, intent(in) :: x_max
  real, intent(inout) :: x

  real :: dx
  real :: s1
  real :: s2
  real :: residual

  real, external :: p32


  residual = 1.0

  dx = x_max/2.0-eps_x
  x  = x_max/2.0

  do while (residual > eps_x)

     s1 = p32(x)
     s2 = p32(x+dx)

     dx = 0.5*dx

     if (s1 <= s .and. s2 >= s) then
        x = x+dx
     else
        x = x-dx
     endif

     residual = abs(s2-s1)

  enddo

end subroutine invert_p32


