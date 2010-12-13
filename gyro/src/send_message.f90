!-----------------------------------------------------------
! send_message.f90
!
! PURPOSE:
!  Routine to write message to runfile.
!
! REVISIONS:
! 31 Jan 06: jc
! 16 Jan 07: mrf
!  "trim" runfile.
!-----------------------------------------------------------

subroutine send_message(message)

  use gyro_globals

  implicit none

  character (len=*), intent(in) :: message

  select case (output_flag)

  case (1)

    if (i_proc == 0) then 
        open(unit=1,file=trim(runfile),status='old',position='append')
        write(1,10) '----------------------------------------------------------'
        write(1,10) message
        write(1,10) '----------------------------------------------------------'
        close(1)
     endif

  end select

10 format(a)

end subroutine send_message
