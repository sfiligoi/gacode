c -*-f90-*-
c
c     This file contains the common block used to control input/output
c     for the ecrh ray tracing code TORAY.
c
         real  fdout,         tout
         integer io1,      io2,           io3,           io4,
     .      io5,           io6,           io7,           io8,
     .      io9,           io0,           iotty,         numout
c
         COMMON /TorGA_comio/
     .      fdout,         tout,
     .      io1,           io2,           io3,           io4,
     .      io5,           io6,           io7,           io8,
     .      io9,           io0,           iotty,         numout,
     .      nlout,         cqldat,        netcdfdat
c
         COMMON /TorGA_comion/
     .      fname1,        fname2,        fname3,        fname4,
     .      fname5,        fname6,        fname7,        fname8,
     .      fname9,        fname0
c
         CHARACTER
     .      fname1*8,      fname2*8,      fname3*8,      fname4*8,
     .      fname5*8,      fname6*8,      fname7*8,      fname8*8,
     .      fname9*8,      fname0*8
c


         LOGICAL
     .      nlout(0:noutmx),   cqldat,  netcdfdat
c
         DIMENSION
     .      tout(noutmx)
c   noutmx is already declared integer in paradm.i
c
