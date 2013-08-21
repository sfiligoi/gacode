      MODULE IONS
      USE param, only : kprim,kimp,kk,kion,kj
! --- old INCLUDE file ions.i
      implicit none
!
!
      logical,public      ::  squeeze
      integer*4,public    ::  adjzeff
      integer*4,public    ::  njenp(kprim),njeni(kimp)
      character*8 ,public ::  namep(kprim)
      character*8 ,public ::  namei(kimp)
      character*8 ,public ::  namen(2)
      character*8 ,public ::  nameb
      character*8 ,public ::  nameu(kk)
      character*8 ,public ::  cer_ion
!
      real *8,public      :: atw(kion)
      real *8,public      :: atomno(kion)
      real *8,public      :: dzdtim(kj, kimp)
      real *8,public      :: z(kj,kion)
      real *8,public      :: zsq(kj,kion)
      real *8,public      :: dzdte(kj, kion)
      real *8,public      :: zeff(kj)
      real *8,public      :: rfatmwt,fd_thermal
!
      END MODULE IONS
