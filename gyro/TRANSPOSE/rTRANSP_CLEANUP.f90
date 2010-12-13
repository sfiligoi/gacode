!---------------------------------------------
! rTRANSP_CLEANUP.f90
!
! PURPOSE:
!  Deallocation of arrays used in rTRANSP
!
! NOTES:
!  Communication routine originally used 
!  in old GYRO.
!
! REVISIONS
! 03 June 02: jc
!  Documented.
!---------------------------------------------

subroutine rTRANSP_CLEANUP

  use rTRANSP_GLOBALS

  implicit none

  deallocate(s) 

  deallocate(ij)
  deallocate(i_ij)
  deallocate(j_ij)

  deallocate(i_rc)
  deallocate(ki)
  deallocate(k_ki)
  deallocate(i_ki)

  deallocate(j_map)
  deallocate(p_ki_loc_map)

  deallocate(q_send)
  deallocate(q_recv)

  return
  
end subroutine rTRANSP_CLEANUP
