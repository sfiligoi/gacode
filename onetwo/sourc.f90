      MODULE sourc
      USE param, only : kk,kj,kion,kprim,krf
      real *8                                                             &
                  siadd(kj,2), s(kk,kj), sione(kj), sion(kj,2),           &
                   sfus(kj), stfus(kj), sbfus(kj), srecom(kj,2),          &
                   scx(kj,2), sbeam(kj), qdimpl(kj), qdelt(kj),           &
                    qgam(kj), qohm(kj), qione(kj), qioni(kj), qcx(kj),    &
                    qbeame(kj), qbeami(kj), qrfe(kj), qrfi(kj),           &
                    qrad(kj), qfus(kj), qfuse(kj), qfusi(kj),             &
                    qbfus(kj), qbfuse(kj), qbfusi(kj), qtfus(kj),         &
                    qtfuse(kj), qtfusi(kj), qmag(kj), qsawe(kj),          &
                    qsawi(kj), ssaw(kj, kprim), ssawe(kj), qexch(kj),     &
                    wdelt, wgam, wohm, refrad, s2d(kj, kion),             &
                    sum2d(kj), qe2d(kj), qi2d(kj), flxmod(2),             &
                    tnion0(2), snadd(2), totohm, totbe, totbi,            &
                    snaddt(2), sn2d, sngas(2), flxadd(2), totrf,          &
                    curohm(kj), curbe(kj), curbi(kj), currf(kj),          &
                    curdri(kj), eta(kj), etap(kj), esaw(kj),              &
                    njqin(krf), qine(kj, krf), qini(kj, krf),             &
                    wejcm3(kj), wijcm3(kj),dpedtc(kj),dpidtc(kj),         &
                    wbjcm3(kj), pohm1(kj), pohm2(kj), cur2(kj),           &
                    conde(kj), conve(kj), condi(kj), convi(kj),           &
                    qconde(kj), qconve(kj), qcondi(kj), qconvi(kj),       &
                    totboot, curboot(kj),curboot_bt(kj),curboot_bp(kj),   &
                    qohmi(kj), curbet(kj), curb_external(kj),             &
                    totb, curb(kj), scurdri(kj), scurfast(kj),            &
                    pconde(kj), pcondi(kj), pconve(kj), pconvi(kj),       &
                    xjbni(kj, kion), xjbti(kj, kion), xjbne(kj),          &
                    xjbte(kj), xjbnf(kj), xlica, gconde(kj),              &
                    gcondi(kj),curdbeam(kj), dudtsv(kk,kj),               &
                    totdri, totbeam,gconve(kj),gconvi(kj),                &
                    gconvde(kj),gconvdi(kj),sbfsav(kj),qbfsav(kj),        &
                    qbeame_intg, qbeami_intg,curb_intg,qbth_intg,         &
                    qbth(kj),wbeam_intg,enbeam_intg,sbeam_intg,           &
                    jtor_eccd(kj,krf),mult_curboot,spellet(kj)

!   spellet is source of elect/ions due to Paul Parks model. Has 
!   nothing to do with propel pellet model at this time.
      Integer,private :: k,j
      INTEGER,PRIVATE,PARAMETER :: ktot = kj*krf
                    Data ((jtor_eccd(j,k),j=1,kj),k=1,krf) /ktot*0.0/
!

      END MODULE sourc
