!-----------------------------------------------------------
! tgyro_multi_driver.f90
!
! PURPOSE:
!  Main driver for multi-job utility.
!----------------------------------------------------------

subroutine tgyro_multi_driver

  use mpi
  use tgyro_globals
  use gyro_interface

  implicit none

  integer :: i

  ! GYRO ERROR STATUS:
  integer, dimension(n_inst) :: gyro_error_status_vec
  character (len=80), dimension(n_inst) :: gyro_error_message_vec

  gyro_restart_method = 1

  ! See gyro/src/gyro_globals.f90 for definition of transport_method
  transport_method = 1

  ! Initialize GYRO
  call gyro_init(paths(color+1),gyro_comm)

  ! Run GYRO
  call gyro_run(gyrotest_flag,gyro_restart_method,transport_method)

  call MPI_ALLGATHER(gyro_error_status_out,1,MPI_INTEGER,&
       gyro_error_status_vec,1,MPI_INTEGER,gyro_adj,ierr)
  call MPI_ALLGATHER(gyro_error_message_out,80,MPI_CHARACTER,&
       gyro_error_message_vec,80,MPI_CHARACTER,gyro_adj,ierr)

  if (i_proc_global == 0) then

     open(unit=1,file=trim(runfile),position='append')
     do i=1,n_inst
        write(1,*) trim(paths(i)),": ",trim(gyro_error_message_vec(i))
     enddo
     close(1)

  endif

end subroutine tgyro_multi_driver

