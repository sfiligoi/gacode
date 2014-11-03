!-----------------------------------------------------------
! cgyro_run.f90
!
! PURPOSE:
!  Manage call to local GKCOLL simulation.
!---------------------------------------------------------

subroutine cgyro_run()

   use cgyro_globals
   use cgyro_interface

   implicit none

   integer :: is

   ! Map INTERFACE parameters -> GLOBAL variables
   call map_interface2global

   ! Run GKCOLL
   call cgyro_do

   cgyro_error_status_out  = error_status
   cgyro_error_message_out = error_message

 end subroutine cgyro_run
