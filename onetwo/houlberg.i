c
c --- INCLUDE file houlberg.i
c
c **********************************************************************
c                                                                      *
c   Common blocks for NCLASS Routine implemented by Jon Kinsey         *
c   and original bootstrap routine implemented by Daniel Finkenthal.   *
c                                                                      *
c   NOTE: The include file 'param.i' used by the rest of ONETWO        *
c         must be included BEFORE this include file in the calling     *
c         program - kprim and kimp set in param.i                      *
c                                                                      *
c   Date Created:       6/4/95                                         *
c   Last Modified:      9/25/99                                        *
c                                                                      *
c **********************************************************************
c
c *** mxns - maximum number of ion      species (Z<3)
c *** mxmh - maximum number of impurity species (Z>2)
c *** mxmi - maximum number of total    species, including electrons
c
      integer         mxns,         mxns2,
     .                mxmh,         mxmi
c
****  parameter      (mxns = kprim, mxns2 = mxns + 2,
**** .                mxmh = kimp , mxmi  = 1 + mxns + mxmh)
c
c *** Add +2 to mxns for beam, alpha ions
c
      parameter      (mxns = kprim + 2, mxns2 =     mxns + 2,
     .                mxmh = kimp     , mxmi  = 1 + mxns + mxmh)
c
c ----------------------------------------------------------------------
c
c *** mxmz - maximum number of ion charge states - transport
c
      integer         mxmz, mx_mi, mx_mz, mx_ms
c      parameter      (mxmz = 10)
       parameter       (mxmz = 8)   ! JEK 10/13/99
       parameter      (mx_mi=mxmi,mx_mz=mxmz,mx_ms=mxns)
c
c ----------------------------------------------------------------------
c
c *** comcon - physics and machine constants
c
      real*8         clight,          coulom,          elemas,
     .               epsilo,          pi,              promas,
     .               xj7kv,           xmachl,          xmachp,
     .               xmachs,          xmuo,            zero
      integer        imachn,          imachp,          lchar,
     .               lword
      common /comcon/clight,          coulom,          elemas,
     .               epsilo,          pi,              promas,
     .               xj7kv,           xmachl,          xmachp,
     .               xmachs,          xmuo,            zero,
     .               imachn,          imachp,          lchar,
     .               lword
c
c ----------------------------------------------------------------------
c
c *** comncl - terms for multispecies neoclassical transport
c
      real*8         amuai,           tai,             vtai,
     .               amntau,          denai,
     .               xzi,             capm,            capn,
     .               calm,            caln,
     .               tauab,
     .               xnud,            xnut,            xnuti,
     .               fex,             uai,
     .               gfl,             qfl,             sqz,
     .               xkban,           xkps,            xmu,
     .               dencut,          xgrp,            xgrt,
     .               xft,             xfhat,
     .               xngrth,          xb2,             xbm2,
     .               xgrbm2,          xedotb,          xjdotb,
     .               etap,            fm
      integer        kboot,           mi,              mz,
     .               kncord
      common /comncl/amuai(mxmi),     tai(mxmi),       vtai(mxmi),
     .               amntau(mxmi,mxmi),
     .               denai(mxmi,mxmz), 
     .               xzi(mxmi,mxmz), capm(3,3,mxmi,mxmi),
     .               capn(3,3,mxmi,mxmi),
     .               calm(3,3,mxmi),  caln(3,3,mxmi,mxmi),
     .               tauab(mxmi,mxmi,mxmz,mxmz),
     .               xnud(mxmi,mxmz), xnut(mxmi,mxmz),
     .               xnuti(3,mxmi,mxmz), fex(3,mxmi,mxmz),
     .               uai(3,2,mxmi,mxmz),
     .               gfl(4,mxmi,mxmz),qfl(4,mxmi,mxmz),
     .               xkban(mxmi,mxmz),xkps(mxmi,mxmz),
     .               xmu(3,3,mxmi,mxmz),      sqz(mxmi,mxmz),
     .               dencut,          xgrp(mxmi,mxmz), xgrt(mxmi),
     .               xft,             xfhat,
     .               xngrth,          xb2,             xbm2,
     .               xgrbm2,          xedotb,          xjdotb,
     .               etap,            fm(3),
     .               kboot,           mi,              mz,
     .               kncord
c
