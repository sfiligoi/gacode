program neo

  use mpi
  use neo_globals

  implicit none

  integer :: ierr
  integer, external :: omp_get_max_threads

  !----------------------------------------------------------------
  ! Query OpenMP for threads
  !
  n_omp = omp_get_max_threads()

  call MPI_INIT(ierr)

  ! Path is cwd:
  path= './'

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  NEO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(NEO_COMM_WORLD,i_proc,ierr)
  call MPI_COMM_SIZE(NEO_COMM_WORLD,n_proc,ierr)
  !-----------------------------------------------------------------

  call neo_read_input
  call neo_do

  call MPI_FINALIZE(ierr)

end program neo
