!--------------------------------------------------------------
! gkcoll_init_serial.f90
!
! PURPOSE:
!  Initialize external GKCOLL interface for serial use.
!---------------------------------------------------------------

subroutine gkcoll_init_serial(path_in)

  use gkcoll_globals
  use gkcoll_interface

  implicit none

  integer :: ierr

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in

  ! Set appropriate global variables
  path = path_in
  i_proc = 0
  n_proc = 1 

  ! This is not a real communicator.
  GKCOLL_COMM_WORLD = -1

end subroutine gkcoll_init_serial
