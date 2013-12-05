!--------------------------------------------------------
! sivukhin.f90
!
! PURPOSE:
!  Compute a low-accuracy but fast approximation to the 
!  ion-alpha heating fraction.
!
!               x
!            1  /     dy
!    F(x) = --- | -----------
!            x  /  1+y^(3/2)
!               0
!          
!  Here, F is the fraction of the alpha energy transferred
!  to ions (at common temperature Ti) by collisions, and
!
!               x = E_alpha/E_crit  
!
!  Details are given in Stix, Plasma Phys. 14 (1972) 367.
!  The function F is derived from Sivukhin's energy loss
!  equation and so that is the rationale for the name.
!---------------------------------------------------------

real function sivukhin(x)

  use math_constants

  implicit none

  real, intent(in) :: x
  integer :: i
  real :: f,dy,yi

  integer, parameter :: n=12


  if (x > 0.1) then

     if (x > 4.0) then

        ! Large-x asymptotic formula

        f = (2*pi/3)/sin(2*pi/3)-2.0/sqrt(x)+0.5/(x*x)

     else

        ! Numerical integration

        dy = x/(n-1)
        f  = 0.0
        do i=1,n
           yi = (i-1)*dy
           if (i == 1 .or. i == n) then
              f = f+0.5/(1.0+yi**1.5) 
           else
              f = f+1.0/(1.0+yi**1.5) 
           endif
        enddo
        f = f*dy

     endif

     sivukhin = f/x

  else

     ! Small-x asymptotic series

     sivukhin = 1.0-0.4*x**1.5

  endif

end function sivukhin
