!-----------------------------------------------------------------
! cgyro_step.f90
!
! PURPOSE:
!  Common work variables for adaptive RK routines
!-----------------------------------------------------------------

module cgyro_step
  
  integer :: itrk,conv

  real :: tol
  real :: orig_delta_t
  real :: deltah2
  real :: deltah2_max, deltah2_min

  real, dimension(2) :: error_sum,error_x
  
  integer, parameter :: itrk_max = 1000
  real, parameter :: eps = 2e-12
  
end module cgyro_step

  
