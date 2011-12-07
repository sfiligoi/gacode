!--------------------------------------------------------------
! gkcoll_init.f90
!
! PURPOSE:
!  Initialize external GKCOLL interface for parallel use.
!---------------------------------------------------------------

subroutine gkcoll_init(path_in,mpi_comm_in)

  use mpi
  use gkcoll_globals
  use gkcoll_interface

  implicit none

  integer :: ierr

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer, intent(in) :: mpi_comm_in

  ! Set appropriate global variables
  path = path_in
  GKCOLL_COMM_WORLD = mpi_comm_in
  call MPI_COMM_RANK(GKCOLL_COMM_WORLD,i_proc,ierr)
  call MPI_COMM_SIZE(GKCOLL_COMM_WORLD,n_proc,ierr)

end subroutine gkcoll_init
