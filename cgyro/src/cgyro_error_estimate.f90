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

  real :: norm
  real :: norm_loc
  real :: error_loc

  if (1 == 1) then

     ! Estimate of field via quadratic interpolation
     field_loc = 3.0*field_old-3.0*field_old2+field_old3

     ! Define norm and error for each mode number n
     norm_loc  = sum(abs(field))
     error_loc = sum(abs(field-field_loc))

     ! Get sum of all errors
     call MPI_ALLREDUCE(error_loc, &
          field_error, &
          1, &
          MPI_DOUBLE_PRECISION, &
          MPI_SUM, &
          NEW_COMM_2, &
          i_err)

     ! Get sum of all norms
     call MPI_ALLREDUCE(norm_loc, &
          norm, &
          1, &
          MPI_DOUBLE_PRECISION, &
          MPI_SUM, &
          NEW_COMM_2, &
          i_err)

  else

     norm_loc  = sum(abs(h_x))
     error_loc = sum(abs(rhs(:,:,1)))

     ! Get sum of all errors
     call MPI_ALLREDUCE(error_loc, &
          field_error, &
          1, &
          MPI_DOUBLE_PRECISION, &
          MPI_SUM, &
          CGYRO_COMM_WORLD, &
          i_err)

     ! Get sum of all norms
     call MPI_ALLREDUCE(norm_loc, &
          norm, &
          1, &
          MPI_DOUBLE_PRECISION, &
          MPI_SUM, &
          CGYRO_COMM_WORLD, &
          i_err)

  endif


  field_error = field_error/norm

  if (error_tol > 0.0) then
     if (field_error > error_tol .and. i_time > 2) then
        call cgyro_error('Integration error exceeded limit.')
     endif
  endif

  field_old3 = field_old2
  field_old2 = field_old
  field_old  = field

end subroutine cgyro_error_estimate
