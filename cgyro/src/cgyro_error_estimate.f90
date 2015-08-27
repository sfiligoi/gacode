!==================================================================================
! Provide integration error estimate via quadratic interpolation.
!==================================================================================

subroutine cgyro_error_estimate

  use cgyro_globals
  use cgyro_io

  implicit none


  ! Estimate of field via quadratic interpolation
  field_loc   = 3.0*field_old-3.0*field_old2+field_old3
  field_error = sum(abs(field-field_loc))/sum(abs(field))

  field_old3 = field_old2
  field_old2 = field_old
  field_old  = field

  if (field_error > 0.1 .and. i_time > 2) then
     call cgyro_error('Integration error > 0.1')
  endif

end subroutine cgyro_error_estimate
