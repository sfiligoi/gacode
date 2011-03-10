!-----------------------------------------------------------------
! gyro.f90
!
! PURPOSE:
!  Main program wrapper. 
!-----------------------------------------------------------------

program gyro

  use gyro_globals

  implicit none

  integer :: ierr

  include 'mpif.h'

  interface
     subroutine gyro_do(skipinit)
       integer, optional :: skipinit
     end subroutine gyro_do
  end interface

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
  transport_method = 0
  !
  call gyro_read_input
  call gyro_read_input_extra
  !
  !----------------------------------------------------------------
  ! Split the communicator if using GKEIGEN added parallelization.
  call GKEIGEN_split_comm
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Reset process ranks if MPI_COMM_WORLD was split.
  call MPI_COMM_RANK(GYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(GYRO_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------

  ! Run gyro.
  call gyro_do
  call send_line('STATUS: '//gyro_exit_message)

  call MPI_FINALIZE(ierr)

end program gyro
