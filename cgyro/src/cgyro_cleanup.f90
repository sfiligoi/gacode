!-----------------------------------------------------------
! gyro_cleanup.f90
!
! PURPOSE:
!  Deallocate and clean up.
!-----------------------------------------------------------

subroutine cgyro_cleanup

  use mpi
  use parallel_lib
  use timer_lib
  use cgyro_globals
  use cgyro_interface

  implicit none

  ! No need for cleanup in test mode
  if (test_flag == 1) return


  call MPI_COMM_FREE(NEW_COMM_1,i_err)
  call MPI_COMM_FREE(NEW_COMM_2,i_err)

  call parallel_lib_cleanup
  call timer_lib_cleanup

  signal = 0
  
end subroutine cgyro_cleanup
