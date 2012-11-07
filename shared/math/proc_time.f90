subroutine proc_time(seconds)

  use mpi
  implicit none

  real :: seconds
 
  seconds = MPI_WTIME()
 
end subroutine proc_time
