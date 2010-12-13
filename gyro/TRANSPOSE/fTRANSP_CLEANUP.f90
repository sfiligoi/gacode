!---------------------------------------------
! fTRANSP_CLEANUP.f90
!
! PURPOSE:
!  Deallocation of arrays used in fTRANSP
!
! NOTES:
!  Communication routine originally used 
!  in old GYRO.
!
! REVISIONS
! 03 June 02: jc
!  Documented.
!---------------------------------------------

subroutine fTRANSP_CLEANUP

  use fTRANSP_GLOBALS

  implicit none

  deallocate(s) 

  deallocate(ij)
  deallocate(i_ij)
  deallocate(j_ij)

  deallocate(i_rc)
  deallocate(jk)
  deallocate(j_jk)
  deallocate(k_jk)

  deallocate(i_map)
  deallocate(p_jk_loc_map)

  deallocate(q_send)
  deallocate(q_recv)

  return
  
end subroutine fTRANSP_CLEANUP
