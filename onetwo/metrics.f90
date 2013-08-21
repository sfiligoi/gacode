   MODULE metrics
      USE mhdpar,only : kpsi
!     Used for Houlberg Bootstrap Model
      IMPLICIT NONE
      REAL *8  gfm(3,kpsi), grth(kpsi),  bsq(kpsi), &
               bmsq(kpsi), grbmsq(kpsi) 
      INTEGER  jhirsh0

   END MODULE metrics
