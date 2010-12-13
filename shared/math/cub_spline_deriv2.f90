!---------------------------------------------------------------
! cub_spline_deriv.f90
!
! PURPOSE:
!  Take known points x(1:n),y(1:n) and perform cubic spline
!  interpolation at the node points xi(1:ni) to give yi(1:ni).
!
!  Use 'natural' cubic spline interpolation, where natural 
!  means the second derivative is zero at the boundaries.
!
!  INPUT  : x(1:n),y(1:n),n
!
!  OUTPUT : yp(1:ni)
!----------------------------------------------------------------

subroutine cub_spline_deriv2(x,y,n,yp)

  !-------------------------------------------------------------
  implicit none
  !
  integer :: i,ii
  real :: x0 
  !
  integer, intent(in) :: n
  real, intent(in), dimension(n) :: x,y
  real, intent(inout), dimension(n) :: yp
  !
  ! LAPACK working variables
  !
  integer :: info
  integer, dimension(n) :: ipiv 
  real, dimension(n) :: c,z
  real, dimension(n-1) :: zl,zu,h
  real, dimension(n-2) :: zu2
  real, dimension(n-1) :: b,d
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Define coefficients of spline matrix 
  !
  do i=1,n-1 
     h(i)  = x(i+1)-x(i)
     zl(i) = h(i)
     zu(i) = h(i)
  enddo
  zl(n-1) = 0.0
  zu(1)   = 0.0

  z(1) = 2*h(1)
  do i=2,n-1
     z(i) = 2.0*(h(i-1)+h(i))
     c(i) = 3.0*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))
  enddo
  z(n) = 2*h(n)
  c(1) = 2*c(2)-c(3)
  c(n) = 2*c(n-1)-c(n-2)
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Solve the system using LAPACK
  !
  call DGTTRF(n,zl,z,zu,zu2,ipiv,info)
  call DGTTRS('N',n,1,zl,z,zu,zu2,ipiv,c,n,info)
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Find remaining polynomial coefficients:
  !
!  c(n) = z(n)
  do i=1,n-1
     b(i) = (y(i+1)-y(i))/h(i)-h(i)*(2.0*c(i)+c(i+1))/3.0
     d(i) = (c(i+1)-c(i))/(3.0*h(i))
  enddo
  !-------------------------------------------------------------

  !print *,c(:)
  
  !-------------------------------------------------------------
  ! Using known polynomial coefficients, perform interpolation.
  !
  !  S(x) = y(i) + b(i) [x-x(i)] + c(i) [x-x(i)]^2 
  !                                        + d(i) [x-x(i)]^3
  !
  ! and then compute derivative based on S'(x). 
  i  = 1
  ii = 1
  do while (ii <= n)
     x0 = x(ii)
     if (x0 <= x(i+1)) then
        yp(ii) = b(i)+(x0-x(i))*(2*c(i)+(x0-x(i))*3*d(i))
        ii = ii+1
     else
        i = i+1
     endif
  enddo
  !-------------------------------------------------------------

end subroutine cub_spline_deriv2

