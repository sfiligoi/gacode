!-----------------------------------------------------------
! tgyro_multi_driver.f90
!
! PURPOSE:
!  Main driver for multi-job utility.
!----------------------------------------------------------

subroutine tgyro_multi_driver

  use tgyro_globals

  implicit none

  integer :: i
  character(len=80), dimension(n_inst) :: ccollect

  include 'mpif.h'

  allocate(gyro_exit_status(n_inst))
  allocate(gyro_exit_message(n_inst))

  !---------------------------------------
  ! Error monitoring variables
  !
  gyro_exit_message(:) = 'INFO: GYRO did not complete'
  !---------------------------------------

  gyro_restart_method = 1
  transport_method    = 1

  ! Initialize GYRO
  call gyro_init(paths(color+1), &
       gyro_comm)

  ! Run GYRO
  call gyro_run(gyrotest_flag, &
       gyro_restart_method, &
       transport_method, &
       gyro_exit_status(color+1), &
       gyro_exit_message(color+1))

  ! Synchronize exit message

  call MPI_ALLGATHER(gyro_exit_message(color+1),&
       80,&
       MPI_CHARACTER,&
       ccollect,&
       80,&
       MPI_CHARACTER,&
       gyro_adj,&
       ierr)

  gyro_exit_message(:) = ccollect(:)

  if (i_proc_global == 0) then

     open(unit=1,file=trim(runfile),position='append')
     do i=1,n_inst
        write(1,*) trim(paths(i)),": ",trim(gyro_exit_message(i))
     enddo
     close(1)

  endif

end subroutine tgyro_multi_driver

