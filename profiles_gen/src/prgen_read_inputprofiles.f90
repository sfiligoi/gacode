!--------------------------------------------------------------
! prgen_read_inputprofiles.f90
!
! PURPOSE:
!  Read input.gacode
!--------------------------------------------------------------

subroutine prgen_read_inputprofiles

  use prgen_globals
  use expro

  implicit none

  expro_ctrl_quasineutral_flag = 0
  expro_ctrl_numeq_flag = 0 
  
  call expro_read('input.gacode')
  
  nx = expro_n_exp

  call allocate_internals

  ! Needed for diagnostic printing
  rmin(:) = expro_rmin(:)
  rmaj(:) = expro_rmaj(:)
       
end subroutine prgen_read_inputprofiles
