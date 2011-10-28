!------------------------------------------------------------
! gyro_write_timers.f90
!
! PURPOSE:
!  Control final calculation and output of code timing data.
!------------------------------------------------------------
 
subroutine gyro_write_timers(datafile,io)

  use mpi
  use gyro_globals

  !-----------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !-------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Initial open

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='replace')
        close(io)

     endif

  case(2)

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='old',position='append')
        close(io)
     endif

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        endfile(io)
        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_timers called]'

end subroutine gyro_write_timers
