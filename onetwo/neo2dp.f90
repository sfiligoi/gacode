       MODULE neo2dp
! --- INCLUDE file neo2dp.i
!
! --- flux surface-related quantities found by subroutine FLUXAV
!
      USE mhdpar,ONLY : kpsi
      IMPLICIT NONE
      REAL*8                                                      &
           epsp(kpsi), xhm2p(kpsi), xi11p(kpsi), xi33p(kpsi),     &
                      xipsp(kpsi), ftnclp(kpsi)
!
    END MODULE neo2dp
