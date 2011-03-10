!-----------------------------------------------------------------
! gyro.f90
!
! PURPOSE:
!  Main program wrapper. 
!-----------------------------------------------------------------

program gyro

  use mpi
  use gyro_globals

  implicit none

  integer :: ierr

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

  !----------------------------------------------------------------
  ! Split the communicator into gkeigen_proc_mult parallel instances
  ! of gyro each calculating a fraction of the columns of the
  ! gkeigen matrix.
  !
  ! j_proc_tot : processor rank in MPI_COMM_WORLD
  ! i_proc     : processor rank in GYRO_COMM_WORLD
  ! n_proc     : (total process)/gkeigen_proc_mult
  !               number of processes in each GYRO_COMM_WORLD
  !
  ! If gkeigen_proc_mult > 1, the user must specify a processor 
  ! number that is an integer N such that:
  !
  !        N = (n_proc)*(gkeigen_proc_mult)
  ! and
  !
  !        n_nek = (n_proc)*M, where M is any integer.
  !
  If (linsolve_method /= 2) gkeigen_proc_mult = 1
  gkeigen_j_set = Mod(i_proc,gkeigen_proc_mult)
  call MPI_COMM_SPLIT(GYRO_COMM_WORLD,&
       gkeigen_j_set,&
       i_proc,&
       GKEIGEN_J_SUBSET, &
       i_err)
  !
  GYRO_COMM_WORLD = GKEIGEN_J_SUBSET
  !
  call MPI_COMM_RANK(GYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(GYRO_COMM_WORLD,n_proc,i_err)

  call MPI_COMM_RANK(MPI_COMM_WORLD,j_proc_tot,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc_tot,i_err)
  !
  !---------------------------------------------------------------

  ! Run gyro.
  call gyro_do
  call send_line('STATUS: '//gyro_exit_message)

  call MPI_FINALIZE(ierr)

end program gyro
