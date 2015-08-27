program cgyro

  use mpi
  use cgyro_globals, only : path, CGYRO_COMM_WORLD, i_proc, n_proc, i_err

  implicit none

  call MPI_INIT(i_err)
  
  ! Path is cwd:
  path= './'

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  CGYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,i_err)
  !-----------------------------------------------------------------

  call cgyro_read_input
  call cgyro_kernel

  call MPI_FINALIZE(i_err)

end program cgyro
