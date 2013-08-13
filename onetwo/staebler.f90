      MODULE staebler

      USE param, ONLY : kj
      IMPLICIT NONE
! --- INCLUDE file staebler.i
!
      REAL *8                                                   & 
                 fs(18), coefa(kj-1), coefb(kj-1), coefc(kj-1), &
                        coefd(kj-1), sperp(kj-1), tstaebler,    &
                        term1(kj-1), term2(kj-1), term3(kj-1)

      CHARACTER(len=8)         staeblrmodl(15)
      INTEGER ifsflag
!


      DATA staeblrmodl                                               &
                        / 'ifsflag', 'aeh'    , 'ael'    , 'aih',    &
                        'ail'    , 'bh'     , 'bl'     , 'ch' ,      &
                       'cl'     , 'alphae' , 'alphai' , 'betah',     &
                        'sigma'  , 'gammah' , 'used ls' /
                                                                


      END MODULE staebler
