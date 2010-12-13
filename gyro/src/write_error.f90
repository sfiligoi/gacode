!-----------------------------------------------------
! write_error.f90
!
! PURPOSE:
!  This routine prints the instantaneous timestep 
!  error.
!
! REVISIONS
! 09 Feb 00: jc
!  Corrections for egyro.
!-----------------------------------------------------

subroutine write_error(datafile,io,action)

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: action
  integer, intent(in) :: io
  character (len=*), intent(in) :: datafile
  !
  integer :: data_loop
  real :: dummy(n_kinetic)
  !---------------------------------------------------


  select case (action)

  case(1)

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(-4,-2,2,4)

     if (i_proc == 0) then

        ! output to screen

        if (silent_flag == 0 .and. linsolve_method == 1) then
           select case (electron_method)

           case (1) 

              print *
              print 10,'t-integration error ({i}): ',time_error(:)

           case (2,4)

              print *
              print 10,'t-integration error ({i},e): ',time_error(:)

           case (3) 

              print *
              print 10,'t-integration error (e): ',time_error(:)

           end select
        endif

        if (output_flag == 1) then

           ! output to file

           open(unit=io,file=datafile,status='old',position='append')
           write(io,20) time_error(:)
           close(io)

        endif

     endif

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')

        do data_loop=0,data_step
           read(io,20) dummy
        enddo

        endfile(io)
        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) print *,'[write_error called]'

10 format(t2,a,6(1pe11.4,1x))
20 format(10(1pe11.4,1x))

end subroutine write_error
