!--------------------------------------------------------------
! neo_init.f90
!
! PURPOSE:
!  Initialize external NEO interface.
!---------------------------------------------------------------

subroutine neo_init(path_in)

   use neo_globals
   use neo_interface

   implicit none

   ! Input parameters (IN) - REQUIRED
   character(len=*), intent(in) :: path_in

   ! Set appropriate global variables
   path = path_in

end subroutine neo_init
