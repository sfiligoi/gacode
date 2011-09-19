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

     if ((i_proc == 0) .and. (gkeigen_j_set==0)) then 
        open(unit=1,file=trim(runfile),status='old',position='append')
        write(1,*) message
        close(1)
     endif

  end select

end subroutine send_line

!-----------------------------------------------------------
! send_message.f90
!
! PURPOSE:
!  Routine to write message to runfile.
!-----------------------------------------------------------

subroutine send_message(message)

  use gyro_globals

  implicit none

  character (len=*), intent(in) :: message

  select case (output_flag)

  case (1)

    if ((i_proc == 0) .and. (gkeigen_j_set == 0)) then 
        open(unit=1,file=trim(runfile),status='old',position='append')
        write(1,10) '----------------------------------------------------------'
        write(1,10) message
        write(1,10) '----------------------------------------------------------'
        close(1)
     endif

  end select

10 format(a)

end subroutine send_message

!-----------------------------------------------------------
! send_message_real.f90
!
! PURPOSE:
!  Routine to write formatted message to runfile.
!-----------------------------------------------------------

subroutine send_message_real(message,x)

  use gyro_globals

  implicit none

  character (len=*), intent(in) :: message
  real, intent(in) :: x

  select case (output_flag)

  case (1)

     if ((i_proc == 0) .and. (gkeigen_j_set == 0)) then 
        open(unit=1,file=trim(runfile),status='old',position='append')
        write(1,10) '----------------------------------------------------------'
        write(1,'(a,f5.3)') message,x
        write(1,10) '----------------------------------------------------------'   
        close(1)
     endif

  end select

10 format(a)

end subroutine send_message_real
