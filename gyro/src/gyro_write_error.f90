!--------------------------------------------------------
! gyro_write_error.f90
!
! PURPOSE:
!  This routine prints the instantaneous timestep error.
!--------------------------------------------------------

subroutine gyro_write_error(datafile,io)

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer :: data_loop
  real :: dummy(n_kinetic)
  !---------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     if (i_proc == 0) then

        ! output to file

        open(unit=io,file=datafile,status='old',position='append')
        write(io,10) time_error(:)
        close(io)

     endif

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')

        do data_loop=0,data_step
           read(io,10) dummy
        enddo

        endfile(io)
        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_error called]'

10 format(10(1pe11.4,1x))

end subroutine gyro_write_error
