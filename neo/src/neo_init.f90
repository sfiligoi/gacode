!--------------------------------------------------------------
! neo_init.f90
!
! PURPOSE:
!  Initialize external NEO interface for parallel use.
!---------------------------------------------------------------

subroutine neo_init(path_in,mpi_comm_in)

  use mpi
  use neo_globals
  use neo_interface

  implicit none

  integer :: ierr

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer, intent(in) :: mpi_comm_in

  ! Set appropriate global variables
  path = path_in
  NEO_COMM_WORLD = mpi_comm_in
  call MPI_COMM_RANK(NEO_COMM_WORLD,i_proc,ierr)
  call MPI_COMM_SIZE(NEO_COMM_WORLD,n_proc,ierr)

end subroutine neo_init
