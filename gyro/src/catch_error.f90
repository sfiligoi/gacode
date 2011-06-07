!-----------------------------------------------------------
! catch_error.f90
!
! PURPOSE:
!  Routine to print error message, finalize MPI, and
!  stop program execution gracefully. 
!-----------------------------------------------------------

subroutine catch_error(message)

  use gyro_globals

  implicit none

  character (len=*), intent(in) :: message

  select case (output_flag)

  case (0)

     if ((i_proc == 0) .and. (gkeigen_j_set == 0)) then
       print *, '----------------------------------------------------------'
       print *, message
       print *, '----------------------------------------------------------'
     endif

  case (1)

     if ((i_proc == 0) .and. (gkeigen_j_set == 0)) then 
        open(unit=1,file=trim(runfile),status='old',position='append')
        write(1,10) '----------------------------------------------------------'
        write(1,10) message
        write(1,10) '----------------------------------------------------------'
        close(1)
     endif

  end select

  call MPI_finalize(i_err)
  stop

10 format(a)

end subroutine catch_error
