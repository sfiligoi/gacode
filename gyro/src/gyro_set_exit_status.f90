!-----------------------------------------------------------
! gyro_set_exit_status.f90
!
! PURPOSE:
!  Set exit status and exit message.
!-----------------------------------------------------------

subroutine gyro_set_exit_status(message,stat)

  use gyro_globals 

  implicit none

  character (len=*), intent(in) :: message
  integer, intent(in) :: stat

  gyro_exit_status  = stat
  gyro_exit_message = trim(message)

end subroutine gyro_set_exit_status
