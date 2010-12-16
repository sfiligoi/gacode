!----------------------------------------------------------
! ORB_lambda.f90
!
! PURPOSE:
!  Simple routine to get values of lambda_tp and lambda_max
!----------------------------------------------------------

subroutine ORB_lambda(lambda_tp,lambda_max)

  use ORB_private

  !-----------------------------------------------------------
  implicit none
  !
  real, intent(inout) :: lambda_tp
  real, intent(inout) :: lambda_max
  !-----------------------------------------------------------

  lambda_max = ORB_lambda_max
  lambda_tp  = ORB_lambda_tp

end subroutine ORB_lambda
