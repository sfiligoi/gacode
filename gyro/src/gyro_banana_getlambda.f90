!----------------------------------------------------------
! gyro_banana_getlambda.f90
!
! PURPOSE:
!  Simple routine to initialize lambda_tp and lambda_max
!----------------------------------------------------------

subroutine gyro_banana_getlambda(lambda_tp_in,lambda_max_in)

  use gyro_banana_private

  !-----------------------------------------------------------
  implicit none
  !
  real, intent(inout) :: lambda_tp_in
  real, intent(inout) :: lambda_max_in
  !-----------------------------------------------------------

  lambda_max_in = lambda_max
  lambda_tp_in  = lambda_tp

end subroutine gyro_banana_getlambda
