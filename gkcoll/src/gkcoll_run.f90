!-----------------------------------------------------------
! gkcoll_run.f90
!
! PURPOSE:
!  Manage call to local GKCOLL simulation.
!---------------------------------------------------------

subroutine gkcoll_run()

   use gkcoll_globals
   use gkcoll_interface

   implicit none

   integer :: is

   ! Map INTERFACE parameters -> GLOBAL variables
   call map_interface2global

   ! Run GKCOLL
   call gkcoll_do

   gkcoll_error_status_out  = error_status
   gkcoll_error_message_out = error_message

 end subroutine gkcoll_run
