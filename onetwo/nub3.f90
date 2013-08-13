     MODULE nub3
     USE param,only : kprim,kz,kimp
     REAL *8          znipm(kprim), atwpm(kprim),                  &
                      atwim(kimp), zniim(kimp), zti(kz) 
!
     INTEGER                                                       &
                      iexcit, ilorent, mstate, ncont, kdene,       &
                      kdeni, kdenz, ksvi, ksvz, ksve, krad, ngh,   &
                      ngl, iz(kimp), izstrp(kimp)
!
! --- nub3 is is used for variables related to hexnb routine
! --- note that nouthx and ncorin, also required for hexnb,
! --- have been added to INCLUDE file io.i
!
!
      END MODULE nub3
