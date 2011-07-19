!--------------------------------------------------------------
! neo_init.f90
!
! PURPOSE:
!  Initialize external NEO interface for serial use.
!---------------------------------------------------------------

subroutine neo_init_serial(path_in)

  use neo_globals
  use neo_interface

  implicit none

  integer :: ierr

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in

  ! Set appropriate global variables
  path = path_in
  i_proc = 0
  n_proc = 1 

  ! This is not a real communicator.
  NEO_COMM_WORLD = -1

end subroutine neo_init_serial
