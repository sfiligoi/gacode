!--------------------------------------------------------------
! EXPRO_read.f90
!
! PURPOSE:
!  Read experimental profiles (serial).
!
! NOTES:
!  Variable naming corresponds to input.profiles 
!  documentation at 
!
!    http://fusion.gat.com/theory/input.profiles
!--------------------------------------------------------------

subroutine EXPRO_read(path)

  use EXPRO_interface

  implicit none

  character(len=*) :: path

  call EXPRO_read_driver(path)

  call EXPRO_compute_derived

end subroutine EXPRO_read

