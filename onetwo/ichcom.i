c
c --- INCLUDE file ichcom.i
c
c --- NOTE that common block /ICHCOM/ is now properly data-aligned.
c               PLEASE maintain this when you make modifications. ... JF
c
      common /ichcom/ xkpar, ykperp, freq, rmrf(kich),
     .                ngrid, nhigh, kinc
c
