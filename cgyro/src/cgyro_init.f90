!--------------------------------------------------------------
! cgyro_init.f90
!
! PURPOSE:
!  Initialize external CGYRO interface for parallel use.
!---------------------------------------------------------------

subroutine cgyro_init(path_in,mpi_comm_in)

  use mpi
  use cgyro_globals
  !use gkcoll_interface

  implicit none

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer, intent(in) :: mpi_comm_in

  ! Set appropriate global variables
  path = path_in
  CGYRO_COMM_WORLD = mpi_comm_in
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,i_err)

end subroutine cgyro_init
