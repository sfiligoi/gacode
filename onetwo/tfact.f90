      Module tfact
      USE param,only : kbctim,kj,ksplin
      INTEGER,parameter :: n_elms = 5
      real *8                                                          &
                    wneo(5,5), wneot, w3cla, wtyp,                     &
                    w1typ(kbctim), w2typ(kbctim),                      &
                    w3typ(kbctim), w4typ(kbctim),                      &
                    vtyp(kbctim), w12typ, w13typ,                      &
                    w1pro(ksplin,kbctim), w1fact(kj),                  &
                    w2pro(ksplin,kbctim), w2fact(kj),                  &
                    w3pro(ksplin,kbctim), w3fact(kj),                  &
                    w4pro(ksplin,kbctim), w4fact(kj),                  &
                    vpro(ksplin,kbctim),  vfact(kj),                   &
                    wsaw, w1saw, w2saw, w3saw, betlim,                 &
                    neomult1,neomult2,qneomult,w1typmin,w2typmin,      &
                    w3typmin,w4typmin,w2cg,w3cg,w2cgm,w3cgm,           &
                    betlim0, w1typt, w2typt, w3typt, w4typt, vtypt,    &
                    relaxtyp, wstd, typa(16), vre(kj), set_chie_chii,  &
                    r_clamp,chi_clamp,r_elm(n_elms),t_elms(n_elms),    &
                    etam_elm(n_elms),t_elme(n_elms), itot_elm(n_elms)
!
     integer                                                           &
                    jgboot, jhirsh, jneo, iboot, ibaloo,               &
                    nw1pro, nw2pro, nw3pro, nw4pro, nvpro,             &
                    neomult,jboot,k_elms
      character*8    ftcalc, resistive
     END MODULE tfact
