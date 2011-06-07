!-----------------------------------------------------------------
! gyro.f90
!
! PURPOSE:
!  Main program wrapper. 
!-----------------------------------------------------------------

program gyro

  use mpi
  use gyro_globals

  !-----------------------------------------------------------------
  implicit none
  !
  integer :: ierr
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(ierr)
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Path is cwd:
  !
  path= './'
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  GYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(GYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(GYRO_COMM_WORLD,n_proc,i_err)
  !
  ! Standard standalone operation
  transport_method = 1
  !
  call gyro_read_input
  call gyro_read_input_extra
  !
  !----------------------------------------------------------------
  ! Split the communicator if using GKEIGEN added parallelization.
  call GKEIGEN_split_comm
  !---------------------------------------------------------------

  ! Run gyro.
  call gyro_do
  call send_line('STATUS: '//gyro_exit_message)

  call MPI_FINALIZE(ierr)

end program gyro
