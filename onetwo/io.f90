      MODULE io
!
!
      USE param, ONLY : krf
! --- INCLUDE file io.i
!
      PARAMETER  (kddebug = 50, nprtlst_max = 350, npltlst_max = 350)
!
      CHARACTER   versid*16, eqdskin*64
      REAL *8                                                             &
                  timprt, timplt,prtlst(nprtlst_max),                     & 
                 pltlst(npltlst_max), xdebug(20),                         &
                  ddebug(kddebug),banktime, extime                      
      INTEGER                                                            &
                 ncrt, nin, nout, nqik, neqplt, ntrplt, nitre, ngreen,   &
                 nbplt, nbdep, neq, nsavsol, nscr, nfw, nrguess,         &
                 nwguess, nmix, ntweak, nupel, nitrex, ndset,            &
                 iprt, mprt,n_rf,nlog_gcnmp,io_temp,                     &
                  jprt, mplot,nb_strt,                                   &
                 jflux, jcoef, jsourc, jbal, jtfus, ihead, ineu, inub,   &
                 irfcalc(krf), ifred,                                    &
                 nterow, jterow(10), ilastp, itimav,                     &
                 nouthx, ncorin, nunadas, nunncl,                        &
                 use_efit_cntr, extend_seval,                            &
                 nmcgo1, nmcgo2, nmcgo3, ntcc, nmhddat,nback,            &
                 n66,n77,save_scratch1,ioftn(10),iopntr(10), nunits,     &
                 lun_nubeam,io_toray_hist,prtlst_used(nprtlst_max,krf)
!
!
      END MODULE io
