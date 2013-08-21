!
      MODULE nub2
      USE param, only : kj,ke,kb,kf,kz,maxp
      implicit none
      integer*4,public    ::  itrapfi, itrapech, ibcur, ibcx,          &
                              iborb, ibslow,inubpat, npat(2), npts
!
      real *8,public,save ::                                            &
                       bencap(kj,ke,kb), fbe(kj,ke,kb), fbi(kj,ke,kb),  &
                     fbth(kj,ke,kb), bke(kj,ke,kb), bki(kj,ke,kb),      &
                     b1ins(kf), b1ots(kf), b2ins(kf), b2ots(kf),        &
                     ecrit(kj), emzrat(kj), enbeam(kj), enbs(kj),       &
                     enb(kj,ke,kb), enbsav(kj,ke,kb), enbav(kj,ke,kb),  &
                     enbav0(kj), enbav1(kj), forb(ke,kb), fb11(ke,kb),  &
                     fb10(ke,kb), fb01(ke,kb), fb00(ke,kb), fber(ke,kb),&
                     bion(ke,kb), bneut(ke,kb), hdep(kj,ke,kb),         &
                     hdepz(kz,ke,kb), pb0(kj,ke,kb),tbeam(kj),          &
                     ppb(kj,ke,kb), ppbsav(kj,ke,kb), ppbav(kj,ke,kb),  &
                     pinsid(kf), potsid(kf), rinsid(kf), rotsid(kf),    &
                     fpsio(kf), qbsav(kj,ke,kb), sbsav(kj,ke,kb),       &
                     fpsii(kf), spbsav(kj,ke,kb), taupb(kj,ke,kb),      &
                     tauppb(kj,ke,kb), taueb(kj,ke,kb), taus(kj),       &
                     wbeam(kj),wb(kj,ke,kb), wbsav(kj,ke,kb),           &
                     wbav(kj,ke,kb), olossc(kj,ke,kb), oloss(kj,ke,kb), &
                     wb11(ke,kb), wb10(ke,kb), wb01(ke,kb), wb00(ke,kb),&
                     xpts(maxp), ypts(maxp), zpts(maxp), rpts(maxp),    &
                     zetaz(kz,ke,kb), wbeamperp(kj),wbeampara(kj),      &
                     nb_thresh ,enbmin ,tmin_curray,enbmin_curray,      &
                     vx(maxp),vy(maxp),vz(maxp),pitch_a(maxp)
      data wbeam /kj*0.0 /
      data enbeam /kj*0.0/
      data enbmin / 1.e3 /                     !min fast ion density
!
      END MODULE nub2
