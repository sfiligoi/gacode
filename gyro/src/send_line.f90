!-----------------------------------------------------------
! send_line.f90
!
! PURPOSE:
!  Routine to write line to runfile.
!-----------------------------------------------------------

subroutine send_line(message)

  use gyro_globals

  implicit none

  character (len=*), intent(in) :: message

  select case (output_flag)

  case (1)

     if (i_proc == 0) then 
        open(unit=1,file=trim(runfile),status='old',position='append')
        write(1,*) message
        close(1)
     endif

  end select

end subroutine send_line
