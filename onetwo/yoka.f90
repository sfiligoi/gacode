      MODULE yoka
      USE param, only :kj
      real *8                                                   &
               poh, pbe, pbi, sthru, pfil, prf, ataue(kj),      &
               tauer(kj), xndinv(kj), xkeinv(kj), xkiinv(kj),   &
               alptix, gamtix, ptor, pradt,                     &
               chietrinv(kj), chiitrinv(kj)                    
      integer                                                   &
               iyoka, ishot, itime
      END MODULE yoka
