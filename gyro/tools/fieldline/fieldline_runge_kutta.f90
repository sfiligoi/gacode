module fieldline_runge_kutta

contains

  ! Second-order Runge-Kutta

  subroutine rk2(y,y0,h,n)

    implicit none

    integer, intent(in) :: n
    real, intent(in) :: h
    real, intent(in), dimension(n) :: y0
    real, intent(inout), dimension(n) :: y

    real, dimension(n) :: yp,k1,k2

    call func(y0,yp,n)
    k1 = h*yp

    call func(y0+k1,yp,n)
    k2 = h*yp

    y = y0+0.5*(k1+k2)

  end subroutine rk2

  ! Fourth-order Runge-Kutta

  subroutine rk4(y,y0,h,n)

    implicit none

    integer, intent(in) :: n
    real, intent(in) :: h
    real, intent(in), dimension(n) :: y0
    real, intent(inout), dimension(n) :: y

    real, dimension(n) :: yp,k1,k2,k3,k4

    call func(y0,yp,n)
    k1 = h*yp

    call func(y0+0.5*k1,yp,n)
    k2 = h*yp

    call func(y0+0.5*k2,yp,n)
    k3 = h*yp

    call func(y0+k3,yp,n)
    k4 = h*yp

    y = y0+(k1+2*k2+2*k3+k4)/6.0

  end subroutine rk4

end module fieldline_runge_kutta
