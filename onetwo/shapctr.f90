
      MODULE shapctr

      USE mhdpar,only: mxtbcmhd

      IMPLICIT NONE
! --- INCLUDE file shapctr.i
!

      character*8       limpos
!
      REAL *8  tolvol, volaray(mxtbcmhd), volwant,                &
               voladj, dvadjmax, volnudge, voladj1, voladj2,      &
               eqtim1, eqtim2, widep, hitep, xlimpos(3), xlim
      INTEGER itvol 
!
      END MODULE shapctr
