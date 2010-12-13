!-----------------------------------------------------
! p32.f90
!
! PURPOSE:
!  Compute the incomplete Gamma (P) function.
!
! NOTES:
!                     x
!              2      /          -t
! p32(x) = --------   | sqrt(t) e   dt == P(3/2,x)
!          sqrt(pi)   /
!                     0
! 
! REVISIONS
! 10 Jan 01: jeff.candy@gat.com
!  Created.
!---------------------------------------------

real function p32(x)

  use math_constants

  implicit none

  real, intent(in) :: x
  real, external :: DERF


  p32 = DERF(sqrt(x))-sqrt(x)*exp(-x)/(0.5*sqrt(pi))

end function p32
