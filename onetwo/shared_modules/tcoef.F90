      MODULE tcoef
#ifdef ONETWO
        USE param,ONLY : kk,kjm1,kion,kj
#else
        USE nf_param,ONLY : kk,kjm1,kion,kj
#endif
        USE nrtype,                                     ONLY : DP,I4B
        IMPLICIT NONE
        INTEGER *4 use_avg_chi,navg_chi
        REAL(DP)                                                          &
                     d(kk,kk,kjm1), dudr(kk,kjm1), dneo(kjm1),            &
                     dbar(kk,kjm1),comp_term(kk),                         &
                    xkeneo(kj), xkineo(kj), xkicla(kjm1),                 &
                    dtyp(kjm1), xketyp(kjm1), xkityp(kjm1),               &
                    dsaw(kjm1), xkesaw(kjm1), xkisaw(kjm1),               &
                    disl(kjm1), xkeisl(kjm1), xkiisl(kjm1),               &
                    xkebal(kjm1), xkibal(kjm1), alphar(kjm1),             &
                    xnuse(kjm1), xnusi(kjm1), xnus(kion, kjm1),           &
                    ftrap(kj), zetaim(kjm1), xkicha(kjm1),                &
                    xkeohk(kjm1), xkepp(kjm1), dinv(kj), xkecar(kjm1),    &
                    xkedom(kjm1), ydebug(kjm1,50), slene(kj),             &
                    slen(kj,kion), slte(kj), slti(kj), slpres(kj),        &
                    shearp(kj), vpinch(kjm1), xkwmatdm(kj, kion),         &
                    etaim(kj,kion), chiwneo(kj), chiwmatm(kj),            &
                    chieinv(kj), chiinv(kj), xkematdm(kj),                &
                    wetai(kion), etaioff, dfion(kjm1), wetaie,            &
                    fimpurty, errchie(kj), errchii(kj), dfast(kjm1),      &
                    xl31n(kion,kjm1), xl31ne(kjm1), xl32n(kion,kjm1),     &
                    xl32ne(kjm1), xden(kjm1), xdet(kjm1),                 &
                    xdin(kion,kjm1), xdit(kion,kjm1), xkangtyp(kjm1),     &
                    xkangtypv(kjm1), xketypv(kjm1), xkitypv(kjm1),        &
                    dtypv(kjm1),xketot(kj),xkitot(kj),xchiitot(kj),       &
                    xchietot(kj),xdchitot(kj),xkangtot(kj),chiineo(kj),   &
                    chii_ms(kj),chie_ms(kj),etaphi_ms(kj),                &
                    chieinv_transp(kj), chiiinv_transp(kj),               &
                    chie_avg(kj),chii_avg(kj),chie_paleo(kj),             &
                    xke_paleo(kj),chid_avg(kj) !jmp.den
!
     END MODULE tcoef
