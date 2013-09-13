      MODULE mhdbcdtn
    

!
! --- boundary conditon quantities for MHD eq. calculations:
! --- flxeqbcd holds psiloop values of psi at the times
! --- given in timeqbcd.
! --- see file mhdpar.i for info regarding nflxeqbcd
! --- psiloop and psincbcd are the interpolated
! --- values from flxeqbcd at the particular time of interest.
! --- psincbcd ===psi-no coil -boundary condition
! --- it is assumed here that kside >nfcoil, otherwise flxeqbcd
! --- may not be dimensioned correctly.
! --- similarly curfcoil holds current value of fcoil currents
! --- which are stored in fcoilcur as a function of time.
! --- the same relationship exists between probeval and expmp2
! --- and ecoilcur, ecurrt;vescur, vescurrt
! --- pcurmhdt=value of pcurmhd at time of interest
! --- similarly for btormhdt, vlopmhdt
! --- the interpolations in time are all done in sub getmhdbc
! --- psilopcl and probeclc are calculated values of
! --- psiloop and probeval respectively.
! --- for the no coils option, psincbcd holds the values of psi on
! --- the edge of the rectangular MHD grid and bdindex is a pointer
! --- used to generate these boundary values.
!

      USE nrtype,                                   ONLY :  I4B,DP
      USE mhdpar,                                   ONLY :  nesum,nfcoil,nsilop,      &
                                                            magpr2,kside,magpri,      &
                                                            mxtbcmhd,nvessel,nstark,  &
                                                            maxobsrv,nflxeqbcd


!        Original size when using greens tables:
         !*** .  probeval(magpri), expmp2(magpri,mxtbcmhd),
         !*** .  psilopcl(nsilop), probeclc(magpri), xchisq,

         !sizes whne not suing greens table:
         REAL(DP)                                                               &
                    psincbcd(kside), psiloop(nsilop),                           &
                    timeqbcd(mxtbcmhd), flxeqbcd(nflxeqbcd,mxtbcmhd),           &
                    curfcoil(nfcoil), fcoilcur(nfcoil,mxtbcmhd),                &
                    probeval(magpr2), expmp2(magpri,mxtbcmhd),                  &
                    ecoilcur(nesum), ecurrt(nesum,mxtbcmhd),                    &
                    vescur(nvessel), vescurrt(nvessel,mxtbcmhd),                &
                    pcurmhd(mxtbcmhd), btormhd(mxtbcmhd), pcurmhdt, btormhdt,   &
                    vloopmhd(mxtbcmhd), vlopmhdt, expfcoil(nfcoil),             &
                    fcoilinv(nfcoil,maxobsrv),deltat_fixed_boundary(mxtbcmhd),  &
                    errpsilp, errmprbe, errfcoil, errecoil, errvescr,           &
                    fwtvescr(nvessel), fwtecoil(nesum),                         &
                    psilopcl(nsilop), probeclc(magpr2), xchisq,                 &
                    fwtpsilp(nsilop), fwtmprbe(magpri), fwtfcoil(nfcoil),       &
                    rvloop, zvloop

        INTEGER(I4b)                                                            &
                     fixfcoil, fitfcur, ifitpsi, ifitprob, iecurr,              &
                     ivessel, itbcmhd, isgngren, mhdmultp, nobsrvd, minchisq,   &
                     bdindex(kside)

        LOGICAL                 use_stark
        INTEGER(I4B)            istark
!
        REAL(DP)                                                                     &
             tstark(nstark,mxtbcmhd), sigstark(nstark,mxtbcmhd),                     &
             rstark(nstark,mxtbcmhd), fwtstark(nstark,mxtbcmhd),                     &
             zstark(nstark,mxtbcmhd), a1stark(nstark,mxtbcmhd),                      &
             a2stark(nstark,mxtbcmhd), a3stark(nstark,mxtbcmhd),                     &
             a4stark(nstark,mxtbcmhd), timestark(mxtbcmhd),                          &
             tstarkexptl(nstark), sigstarkexptl(nstark),                             &
             rstarkexptl(nstark), zstarkexptl(nstark),                               &
             tstarkcalc(nstark), chisqstark(nstark), chisqstarktot

!
      END MODULE mhdbcdtn
