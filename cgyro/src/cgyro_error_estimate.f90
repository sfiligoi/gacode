!------------------------------------------------------------------------------
! cgyro_error_estimate.f90
!
! PURPOSE:
!  Compute two time-integration error estimates based on 
!  (1) difference between computed field and quadratic extrapolation, 
!  (2) 3rd-order estimate of collisionless error based on RK algebra
!------------------------------------------------------------------------------

subroutine cgyro_error_estimate

  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none

  real, dimension(2) :: norm
  real, dimension(2) :: norm_loc
  real, dimension(2) :: error_loc

  ! 1. Estimate of total (field) error via quadratic interpolation

  field_loc = 3.0*field_old-3.0*field_old2+field_old3

  ! Define norm and error for each mode number n
  norm_loc(1)  = sum(abs(field))
  error_loc(1) = sum(abs(field-field_loc))

  ! 2. Estimate of collisionless error via 3rd-order linear estimate

  norm_loc(2)  = sum(abs(h_x))
  error_loc(2) = sum(abs(rhs(:,:,1)))

  ! Get sum of all errors
  call MPI_ALLREDUCE(error_loc, &
       integration_error, &
       2, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  ! Get sum of all norms
  call MPI_ALLREDUCE(norm_loc, &
       norm, &
       2, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  integration_error = integration_error/norm

  if (integration_error(2) > error_tol .and. i_time > 2) then
     call cgyro_error('Integration error exceeded limit.')
  endif

  field_old3 = field_old2
  field_old2 = field_old
  field_old  = field

end subroutine cgyro_error_estimate
