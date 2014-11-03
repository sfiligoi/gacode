!-----------------------------------------------------------
! catch_error.f90
!
! PURPOSE:
!  Routine to print error message, finalize MPI, and
!  stop program execution gracefully. 
!-----------------------------------------------------------

subroutine catch_error(message)

  use cgyro_globals

  implicit none

  character (len=*), intent(in) :: message

  select case (output_flag)

  case (0)

     if (i_proc == 0) then
       print *, '----------------------------------------------------------'
       print *, message
       print *, '----------------------------------------------------------'
     endif

  case (1)

     if (i_proc == 0) then 
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

!------------------------------------------------
! catch_halt_signal.f90
!
! PURPOSE:
!  halt code on signal.
!------------------------------------------------

subroutine catch_halt_signal

  use mpi
  use cgyro_globals

  !--------------------
  implicit none
  !
  integer :: signal
  logical :: lfe
  !--------------------

  !------------------------------------------------
  ! Halt signal
  !
  if (i_proc == 0) then
     inquire(file='halt',exist=lfe)
     if (lfe .eqv. .true.) then
        open(unit=1,file='halt',status='old')
        read(1,*) signal
        close(1)
     else
        signal = 0
     endif
  endif
  !
  call MPI_BCAST(signal,&
       1,&
       MPI_INTEGER,&
       0,&
       CGYRO_COMM_WORLD,&
       i_err)  
  !------------------------------------------------

end subroutine catch_halt_signal
