      SUBROUTINE tglf_shutdown
!
! save output and dealocate memory 
!
      IMPLICIT NONE

      CALL tglf_harvest
      CALL tglf_deallocate
 
      END SUBROUTINE tglf_shutdown
