!------------------------------------------------------------------------------
! cgyro_error_estimate.f90
!
! PURPOSE:
!  Compute crude time-integration error estimate based on difference
!  between computed field and quadratic extrapolation.
!------------------------------------------------------------------------------

subroutine cgyro_error_estimate

  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none

  real :: field_error_loc

  ! Estimate of field via quadratic interpolation
  field_loc       = 3.0*field_old-3.0*field_old2+field_old3
  field_error_loc = sum(abs(field-field_loc))/sum(abs(field))/n_toroidal

  call MPI_ALLREDUCE(field_error_loc, &
       field_error, &
       1, &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)

  field_old3 = field_old2
  field_old2 = field_old
  field_old  = field

  if (field_error > 0.5 .and. i_time > 2) then
     call cgyro_error('Integration error > 0.5')
  endif

end subroutine cgyro_error_estimate
