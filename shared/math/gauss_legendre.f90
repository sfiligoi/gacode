!---------------------------------------------------------
! gauss_legendre.f90
!
! PURPOSE:
!  Given the lower and upper limits of integration, x1 and 
!  x2, and given n, this routine returns arrays x and w of 
!  length n, containing the abscissas and weights of the 
!  Gauss-Legendre n-point quadrature formula.
!
! NOTES:
!  Taken from "Numerical Recipes" routine GAULEG.
!
!  Book uses eps = 3e-14 -- this has been decreased.
!
! REVISIONS
! 06 Dec 00: jeff.candy@gat.com
!  Copied from book.  
!  Recoded in f90 (with *explicit* variable declarations).
!  Added exception for n=1.
!---------------------------------------------------------

subroutine gauss_legendre(x1,x2,x,w,n)

  implicit none

  real, parameter :: eps = 1e-15

  real, parameter :: pi = 3.141592653589793

  real, intent(in) :: x1
  real, intent(in) :: x2
  integer, intent(in) :: n

  real, dimension(n) :: x
  real, dimension(n) :: w

  real :: xm
  real :: xl
  real :: z
  real :: z1
  real :: p1
  real :: p2
  real :: p3
  real :: pp

  integer :: m
  integer :: j
  integer :: i

  xm = 0.5*(x2+x1)
  xl = 0.5*(x2-x1)

  ! Exception for n=1 is required:

  if (n == 1) then
     x(1) = xm 
     w(1) = 2.0*xl 
     return
  endif

  ! Roots are symmetric.  We need only find half of them.

  m = (n+1)/2

  ! Initialize to fail first do test
  
  z1 = -1.0

  ! Loop over roots.

  do i=1,m

     z = cos(pi*(i-0.25)/(n+0.5))

     do while (abs(z-z1) > eps)

        p1 = 1.0
        p2 = 0.0

        do j=1,n
           p3 = p2
           p2 = p1
           p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
        enddo

        ! p1 is the Legendre polynomial.  Now compute its 
        ! derivative, pp.

        pp = n*(z*p1-p2)/(z*z-1.0)
        z1 = z
        z  = z1-p1/pp

     enddo

     x(i)     = xm-xl*z
     x(n+1-i) = xm+xl*z
     w(i)     = 2.0*xl/((1.0-z*z)*pp*pp)
     w(n+1-i) = w(i)

  enddo

end subroutine gauss_legendre


