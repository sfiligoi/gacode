!--------------------------------------------------
! bound_deriv.f90
!
! PURPOSE:
!  A trivial routine to compute a finite-difference 
!  derivative of a function (not necessarily with 
!  equally-spaced gridpoints).
!
!  The main point is to isolate this oft-used 
!  operation.
!
! REVISIONS:
! 04 April 02: jc 
!  Documented.
!--------------------------------------------------

subroutine per_deriv(df,f,r,n)

  implicit none

  integer, intent(in) :: n

  real, intent(inout), dimension(n) :: df
  real, intent(in), dimension(n) :: f
  real, intent(in),dimension(n) :: r
  real :: dr

  integer :: i

  dr = r(3)-r(1)

  df(1) = (f(2)-f(n))/dr

  do i=2,n-1
     df(i) = (f(i+1)-f(i-1))/dr
  enddo

  df(n) = (f(1)-f(n-1))/dr

end subroutine per_deriv
