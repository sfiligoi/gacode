
  MODULE etc
    USE param, ONLY : kj
    USE nrtype, ONLY : DP
    IMPLICIT NONE
!
     REAL(DP), dimension(:) ::  bp(kj), xtemp(kj), ytemp(kj)
     REAL(DP)  timold, tocur
!
  END MODULE etc
