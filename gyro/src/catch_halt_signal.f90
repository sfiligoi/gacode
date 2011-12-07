!------------------------------------------------
! catch_halt_signal.f90
!
! PURPOSE:
!  halt code on signal.
!------------------------------------------------

subroutine catch_halt_signal

  use mpi
  use gyro_globals

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
       GYRO_COMM_WORLD,&
       i_err)  
  !
  if (abs(signal) == 1) then
     call gyro_set_exit_status('user halt',1)
  endif
  !------------------------------------------------

end subroutine catch_halt_signal
