      SUBROUTINE tglf_shutdown
!
! save output and dealocate memory 
!
      IMPLICIT NONE

      CALL tglf_harvest_local

      CALL tglf_deallocate
 
      END SUBROUTINE tglf_shutdown
