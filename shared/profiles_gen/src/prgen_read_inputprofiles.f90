!--------------------------------------------------------------
! prgen_read_inputprofiles.f90
!
! PURPOSE:
!  Read input.profiles.gen
!--------------------------------------------------------------

subroutine prgen_read_inputprofiles

  use prgen_globals
  use EXPRO_interface

  implicit none

  EXPRO_ctrl_density_method = 1
  EXPRO_ctrl_z(1:3) = 1.0 
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_signq = (-1)*(-1)
  EXPRO_ctrl_signb = -(-1)
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_alloc('./',1) 
  call EXPRO_read

  nx = EXPRO_n_exp

  call allocate_internals

end subroutine prgen_read_inputprofiles
