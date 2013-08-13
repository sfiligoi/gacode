c
c --- INCLUDE file pelcom.i
c
c --- NOTE that common block /PELCOM/ is now properly data-aligned.
c               PLEASE maintain this when you make modifications. ... JF
c
      character       nampel*8
c
      common /pelcom/ pelrad, vpel, timpel(10),
     .                ipel, nbgpel, ipelet, npel,
     .                nampel, pelmod
c
