      MODULE NEUT
      USE param, only : kj,kbctim
!
!
!
! --- INCLUDE file neut.i
!
!
      logical,public ::       reqsp
!
     real *8,public :: reflec
     real *8,public :: gasflx(kbctim,2)
     real *8,public :: twall, wion, wrad, recyc(2)
     real *8,public :: widths, enes, tes
     real *8,public :: relneu, enneu(kj,2), tineu(kj), eirate(kj)
     real *8,public :: eionr(kj), cexr(kj,2), cx12r(kj), volsn(kj,2)
     real *8,public :: enn(kj,2), tn(kj,2), f1w, f2w, fvn
     real *8,public :: ennw(kj,2), ennv(kj,2), vneu(kj), raneut
     real *8,public :: englstn(50), rtandn, rhdn
     real *8,public :: spflux(50,2), fluxn(2), flxmin
     real *8,public :: dn1(kj,2), dn2(kj,2), dnv(kj,2)
     real *8,public :: wn1(kj,2), wn2(kj,2), wnv(kj,2), vneut(kj,2)
     integer *4,public ::in, ineut(2), ipcons(2)
     integer *4,public ::ipfail(2), idiagn, nengn
!
!
     END MODULE  NEUT
