!-----------------------------------------------------
! ORB_s2lambda.f90
!
! PURPOSE:
!  Given s, calculate lambda.
! 
! NOTES:
!  Speed not esential, so use bisection for 
!  stability.
!
! REVISIONS
! 04 Dec 01: jc
!  Documented.
!---------------------------------------------

subroutine ORB_s2lambda(s,lambda)

  use ORB_private

  !-----------------------------------------------------------
  implicit none
  !
  real, parameter :: eps_lambda = 1e-9
  real, intent(in) :: s
  real, intent(inout) :: lambda
  real :: dlambda
  real :: s1
  real :: s2
  real :: residual
  !-----------------------------------------------------------

  residual = 1.0

  dlambda = ORB_lambda_max/2.0-eps_lambda
  lambda  = ORB_lambda_max/2.0

  do while (residual > eps_lambda)

     call ORB_s(lambda,s1)
     call ORB_s(lambda+dlambda,s2)

     dlambda = 0.5*dlambda

     if (s1 <= s .and. s2 >= s) then
        lambda = lambda+dlambda
     else
        lambda = lambda-dlambda
     endif

     residual = abs(s2-s1)

  enddo

end subroutine ORB_s2lambda
