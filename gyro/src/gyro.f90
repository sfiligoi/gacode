!-----------------------------------------------------------------
! gyro.f90
!
! PURPOSE:
!  Main program wrapper for *standalone* GYRO usage.
!-----------------------------------------------------------------

program gyro

  use mpi
  use gyro_globals
  use ompdata

  !-----------------------------------------------------------------
  implicit none
  !
  integer :: ierr
  integer :: supported
  integer, external :: omp_get_max_threads
  !-----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Query OpenMP for threads
  !
  n_omp = omp_get_max_threads()
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator, including support for 
  ! funneled threading (needed if OpenMP is enabled).
  !
  if (n_omp > 1) then
     call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,ierr)
     if (supported < MPI_THREAD_FUNNELED) then
        call catch_error('ERROR: (GYRO) Multi-threaded MPI not supported.')
     endif
  else 
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,ierr)
  endif
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Path is cwd:
  !
  path= './'
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Query MPI for dimensions
  !
  GYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(GYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(GYRO_COMM_WORLD,n_proc,i_err)
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Standard standalone operation
  transport_method = 1
  !
  call gyro_read_input
  call gyro_read_input_extra
  !
  !-----------------------------------------------------------------
  ! initialize OpenMP runtime parameters including block indices
  call gyro_init_ompdata()
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
