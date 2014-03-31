      SUBROUTINE tglf_error(itype,error_message)
!
! error handeling routine for TGLF
! simply writes the error or informatin message
! if itype = 1 the error was fatal so TGLF is shutdown
!
      IMPLICIT NONE
      INTEGER itype
      CHARACTER error_message(32)
!      
      write(0,*)error_message
      if(itype.eq.1)call tglf_shutdown
      RETURN
      END SUBROUTINE tglf_error
!
     LOGICAL FUNCTION tglf_isnan(x)
! 
! test real variables for NAN 
!
     IMPLICIT NONE
     REAL x
     LOGICAL tglf_isnan,isnan

     tglf_isnan = .FALSE.
     if(x.ne.x)tglf_isnan = .TRUE.
     
     RETURN
     END FUNCTION tglf_isnan
!
     LOGICAL FUNCTION tglf_isinf(x)
! 
! test real variables for NAN 
!
     IMPLICIT NONE
     REAL x
     LOGICAL tglf_isinf,isinf

     tglf_isinf = .FALSE.
     if(0.0*x.ne.0.0)tglf_isinf = .TRUE.
     
     RETURN
     END FUNCTION tglf_isinf
   
      
