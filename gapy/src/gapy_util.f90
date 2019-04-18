module util

contains

  !-------------------------------------------------------
  ! bound_deriv.f90
  !
  ! PURPOSE:
  !  Compute the finite-difference derivative of a 
  !  function (on a grid which may be unequally-spaced)
  !  using the Lagrange interpolating polynomial.  
  !-------------------------------------------------------

  subroutine bound_deriv(df,f,r,n)

    implicit none

    integer, intent(in) :: n

    double precision, intent(inout), dimension(n) :: df
    double precision, intent(in), dimension(n) :: f
    double precision, intent(in),dimension(n) :: r

    double precision :: r1,r2,r3,ra
    double precision :: f1,f2,f3

    integer :: i

    do i=1,n

       if (i == 1) then

          ! Left boundary

          ra = r(1)

          r1 = r(1)
          r2 = r(2)
          r3 = r(3)
          f1 = f(1)
          f2 = f(2)
          f3 = f(3)

       else if (i == n) then

          ! Right boundary

          ra = r(n)

          r1 = r(n-2)
          r2 = r(n-1)
          r3 = r(n)
          f1 = f(n-2)
          f2 = f(n-1)
          f3 = f(n)

       else

          ! Interior

          ra = r(i)

          r1 = r(i-1)
          r2 = r(i)
          r3 = r(i+1)
          f1 = f(i-1)
          f2 = f(i)
          f3 = f(i+1)

       endif

       ! Derivative of Lagrange interpolating polynomial:

       df(i) = ((ra-r1)+(ra-r2))/(r3-r1)/(r3-r2)*f3 &
            + ((ra-r1)+(ra-r3))/(r2-r1)/(r2-r3)*f2 &
            + ((ra-r2)+(ra-r3))/(r1-r2)/(r1-r3)*f1

    enddo

  end subroutine bound_deriv

  !-------------------------------------------------------
  ! bound_extrap.f90
  !
  ! PURPOSE:
  !  Extrapolate to left and right boundaries. 
  !-------------------------------------------------------

  subroutine bound_extrap(fa,fb,f,r,n)

    implicit none

    integer, intent(in) :: n

    double precision, intent(inout) :: fa
    double precision, intent(inout) :: fb
    double precision, intent(in), dimension(n) :: f
    double precision, intent(in), dimension(n) :: r

    double precision :: r1,r2,r3
    double precision :: ra,rb
    double precision :: f1,f2,f3


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

end module util
