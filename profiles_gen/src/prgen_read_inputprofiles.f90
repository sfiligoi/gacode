!--------------------------------------------------------------
! prgen_read_inputprofiles.f90
!
! PURPOSE:
!  Read input.profiles
!--------------------------------------------------------------

subroutine prgen_read_inputprofiles

  use prgen_globals
  use vpro

  implicit none

  integer :: i

  expro_ctrl_quasineutral_flag = 0
  expro_ctrl_numeq_flag = 0 
  
  call vpro_read('./')
  
  nx    = expro_n_exp
  n_ion = expro_n_ion

  call allocate_internals

  ! Needed for diagnostic printing
  rmin(:) = expro_rmin(:)
  rmaj(:) = expro_rmaj(:)
       
end subroutine prgen_read_inputprofiles
