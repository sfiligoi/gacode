   MODULE monitr


   USE param, ONLY : kj,kk

   IMPLICIT NONE
!
! --- These arrays are only used to monitor the solution
! --- as it progresses.
!
!
! --- imontf     is total convergence failures
! --- imonn      is failures due to particle conservation
! --- rmonu(i,j) is failure at grid point i due to profile j
!
      INTEGER, PARAMETER :: kmon = kj * kk
      INTEGER    imonn, imontf
      REAL*8,DIMENSION(:,:) ::        rmonu (kj,kk) 
      DATA rmonu / kmon*0.0D0 / ;  DATA   imonn, imontf /0, 0/


   END MODULE monitr
