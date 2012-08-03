!------------------------------------------------
! catch_blowup.f90
!
! PURPOSE:
!  This routine traps a blow-up of the numerical
!  solution by monitoring the potential amplitude. 
!------------------------------------------------

subroutine catch_blowup

  use mpi
  use gyro_globals

  !---------------------------------------------------
  implicit none
  !
  real, parameter :: phi_limit = 0.5
  real :: phi_local
  real :: phi_max
  !---------------------------------------------------

  !------------------------------------------------
  ! Halt on large field
  !
  phi_local = maxval(abs(field_blend))
  !
  call MPI_ALLREDUCE(phi_local, &
       phi_max, &
       1, &
       MPI_DOUBLE_PRECISION, &
       MPI_MAX, &
       NEW_COMM_2, &
       i_err)
  
  if (phi_max > phi_limit .and. nonlinear_flag == 1) then
     call gyro_set_exit_status('field blowup',1)
  endif
  !
  !------------------------------------------------

end subroutine catch_blowup

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
