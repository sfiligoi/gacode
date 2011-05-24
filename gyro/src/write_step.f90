!------------------------------------------------
! write_step.f90 [caller gyro_write_master]
!
! PURPOSE:
!  Output of timestep data.
!------------------------------------------------
 
subroutine write_step(datafile,io)

  use gyro_globals

  !------------------------------------------
  implicit none
  !
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer :: data_loop
  integer :: i_dummy
  real :: dummy
  !------------------------------------------

  select case (io_control)

  case(0)

    return

  case(1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2)

     ! Append

     if (silent_flag == 0 .and. linsolve_method == 1) then
        print *
        print 10,'|---------- Step',step,' of', &
             nstep,'----------|   t = ',t_current
     endif

     open(unit=io,file=datafile,status='old',position='append')
     write(io,20) data_step,t_current
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,20) i_dummy,dummy
     enddo

     endfile(io)
     close(io)

  end select

10 format(t2,a,1x,i6,a,i6,a,f8.3)
20 format(i6,1x,es12.5)

end subroutine write_step
