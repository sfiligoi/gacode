program gkcoll

  use mpi
  use gkcoll_globals, only : path, GKCOLL_COMM_WORLD, i_proc, n_proc

  implicit none

  integer :: ierr

  call MPI_INIT(ierr)
  
  ! Path is cwd:
  path= './'

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  GKCOLL_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(GKCOLL_COMM_WORLD,i_proc,ierr)
  call MPI_COMM_SIZE(GKCOLL_COMM_WORLD,n_proc,ierr)
  !-----------------------------------------------------------------

  call gkcoll_read_input
  call gkcoll_do

  call MPI_FINALIZE(ierr)

end program gkcoll
