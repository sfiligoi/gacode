!
      MODULE machin
      USE param,only : kbctim
! --- the machine id, machinei, is kept along with codeid
! --- in file geom.f90
!
      real *8         rmajor, rminor, elong(kbctim), volume, btor,   &
                     flim, rin, rout, rmin, rmax, zmin, zmax, zax,   &
                     kappa, dkapdt, ltest_code
!
      END MODULE machin
