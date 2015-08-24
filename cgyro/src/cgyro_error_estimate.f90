!==================================================================================
! Provide integration error estimate via quadratic interpolation.
!==================================================================================

subroutine cgyro_error_estimate

  use cgyro_globals

  implicit none

  if (i_time == 1) then

     field_old2 = 0.0

  else 

     ! Estimate of field via quadratic interpolation
     field_loc   = 3.0*field_old-3.0*field_old2+field_old3
     field_error = sum(abs(field-field_loc))/sum(abs(field))

  endif

  field_old3 = field_old2
  field_old2 = field_old
  field_old  = field

end subroutine cgyro_error_estimate
