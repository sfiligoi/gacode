 MODULE soln2d
      USE mhdpar,only: mxtbcmhd
      INTEGER,PARAMETER :: keqmax = mxtbcmhd
!
!
      REAL *8        gcomment(8), dteq, dtmine, deqlst(keqmax),     &
                     vloopcl(keqmax), delcap, delrho, derwght,      &
                     cap_mult, toleq, tolcur, omeq, omcur,          &
                     rmagax, zmagax, rscale, zscale,rbsaxis
     INTEGER*4                                                      &
               irguess, iwguess, ifixshap, ieqmax, ieq, icon,       &
               ibypas, igoitr, itre, maxitr,                        &
               mhdonly, iteq, itcur, ieqfail,ieqdsk, ieqprt,        &
               iprtit, j2prt,rf_output, renormalize_rbp
!
      logical         implicit_fh, xi_include
!
 END MODULE soln2d
