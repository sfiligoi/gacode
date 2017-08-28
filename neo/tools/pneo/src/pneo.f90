program pneo

  use mpi
  use pneo_globals
  use neo_interface
  use EXPRO_interface

  implicit none

  integer :: i,j
  integer :: ni,nj

  !---------------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator.
  !
  call MPI_INIT(i_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,i_err)
  !---------------------------------------------------------------------

  ! Path is cwd:
  path= './'

  ! Read vgen control parameters

  call neo_init(path,MPI_COMM_WORLD)
  neo_silent_flag_in = 1

  ni=4
  nj=4

  do i=1,ni
     do j=1,nj

        neo_q_in = 2.0
        print *,i,j
        call neo_run()
        print *, neo_jpar_dke_out
        
     enddo
  enddo

  call MPI_finalize(i_err)

end program pneo
