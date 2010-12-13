!------------------------------------------------------
! write_local_real.f90
!
! PURPOSE:
!  This routine write a vector of nondistributed reals.
!------------------------------------------------------

subroutine write_local_real(datafile,io,action,n_fn,fn)

  use gyro_globals, only : &
       data_step

  !---------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: action
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  real :: dummy(n_fn)
  !---------------------------------------------------

  select case (action)

  case(1)

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2,4)

     open(unit=io,file=datafile,status='old',position='append')
     write(io,'(20(1pe15.8,1x))')  fn(:)
     close(io)

  case(3)

     ! Reposition after restart

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,*) dummy(:)
     enddo

     endfile(io)
     close(io)

  end select

end subroutine write_local_real
