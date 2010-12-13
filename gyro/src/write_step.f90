!------------------------------------------------
! write_step.f90 [caller gyro_write_master]
!
! PURPOSE:
!  Output of timestep data.
!------------------------------------------------
 
subroutine write_step(datafile,io,action)

  use gyro_globals

  !------------------------------------------
  implicit none
  !
  integer, intent(in) :: action
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer :: data_loop
  integer :: i_dummy
  real :: dummy
  !------------------------------------------

  select case (action)

  case(1)

     ! Initial open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(-4,-2,2,4)

     ! Output

     if (silent_flag == 0 .and. linsolve_method == 1) then
        print *
        print 10,'|---------- Step',step,' of', &
             nstep,'----------|   t = ',t_current
     endif

     if (output_flag == 1) then

        open(unit=io,file=datafile,status='old',position='append')
        write(io,20) data_step,t_current
        close(io)

     endif

     ! Timing information

  case(3)

     ! Reposition after restart

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
