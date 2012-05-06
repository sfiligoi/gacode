!-----------------------------------------------------------
! gyro_kinetic_advance.f90
!
! PURPOSE:
!  Manage implicit advance of fields and electron
!  distribution.  Explicit ion RHS done with gyro_rhs_total.
!-----------------------------------------------------------

subroutine gyro_kinetic_advance

  use gyro_globals

  !-------------------------
  implicit none
  !-------------------------

  call gyro_timer_in('Implicit-he')
  call gyro_get_delta_he
  call gyro_timer_out('Implicit-he')

  call gyro_field_solve_implicit
  call gyro_field_interpolation

  call gyro_timer_in('Implicit-he')
  call gyro_get_he_implicit
  call gyro_timer_out('Implicit-he')

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[gyro_kinetic_advance done]'
  endif

end subroutine gyro_kinetic_advance
