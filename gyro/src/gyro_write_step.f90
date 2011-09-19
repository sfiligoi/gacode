!------------------------------------------------
! gyro_write_step.f90
!
! PURPOSE:
!  Output of timestep data.
!------------------------------------------------
 
subroutine gyro_write_step(datafile,io)

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

        if (nonlinear_flag == 0) then
           print '(a,1pe9.3,a,1pe10.3,1pe10.3,a,1pe9.3,a,5(1pe9.3,1x))',&
                '[t = ',t_current,&
                '][w = ',omega_linear(1,1),&
                '][dw = ',abs(omega_linear(1,2)),&
                '] t_err: ',time_error(:)
        else
           print '(a,1pe9.3,a,5(1pe9.3,1x))',&
                '[t = ',t_current,&
                '] t_err: ',time_error(:)
        endif

     endif

     open(unit=io,file=datafile,status='old',position='append')
     write(io,10) data_step,t_current
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,10) i_dummy,dummy
     enddo

     endfile(io)
     close(io)

  end select

10 format(i6,1x,es12.5)

end subroutine gyro_write_step
