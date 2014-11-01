program cgyro

  use mpi
  use cgyro_globals, only : path, CGYRO_COMM_WORLD, i_proc, n_proc

  implicit none

  integer :: ierr

  call MPI_INIT(ierr)
  
  ! Path is cwd:
  path= './'

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  CGYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,ierr)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,ierr)
  !-----------------------------------------------------------------

  call cgyro_read_input
  call cgyro_do

  call MPI_FINALIZE(ierr)

end program cgyro
