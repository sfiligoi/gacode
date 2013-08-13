c
c --- INCLUDE file mhdpar.i (65x65 version)
c
c --- NOTE: PARAMETER statements also appear in FREYA and NEWBOXX
c --- NOTE: kpsi must be greater than or equal to max of nw and kj
c --- NOTE: kside1,nstark1: for fixed boundary cases mxtbcmhd can
c           be set to a large number ( to accomodate a large list
c           of times in timeqbcd) Consequently the experimental
c           data arrays associated with the free boundary calculations
c           become very large. However  if mxtbcmhd is large enough
c           then we start using virtual memory which destroys the
c           performance of the code to an unacceptable degree.
c           If free boundary equilibria are to be done then
c           set keqmax=mxtbcmhd <~ 300 and set nstark = nstark1
c           set nflxeqbcd = kside, and set magpri  = m
c
      integer    nw, nh, nwh, nfcoil, nrstep, kcoilq, kside, nsilop,
     .           nvessel, nrogow, magpr322, magpri67, magpri, magpr2,
     .           maxobsrv, nesum, mxtbcmhd, kpsi, nstark, necoil,
     .           nstark1, nflxeqbcd, magpri1
c
      parameter (kpsi     =  201                     , ! MAX (nw, kj)
     .           nw       =  129                     , ! WIDTH
     .           nh       =  129                     , ! HEIGHT
     .           nwh      = nw * nh                 ,
     .           nfcoil   =  18                     ,
     .           nrstep   =   1                     ,
     .           kcoilq   = nfcoil * nfcoil         ,
     .           kside    =   2 * (nh + nw - 2)     ,
     .           nflxeqbcd = kside                  ,
     .           nstark1  =   8                     ,
     .           nstark   =   nstark1, ! use for fixed boundary cases
     .           nsilop   =  41                     ,
     .           nvessel  =  24                     ,
     .           nrogow   =   1                     ,
     .           magpr322 =  60                     ,
     .           magpri67 =  29                     ,
     .           magpri1  = magpr322 + magpri67     ,
     .           magpri   = magpri1                 ,
     .           magpr2   = magpr322                ,
     .           maxobsrv = nsilop + magpr2 + nfcoil,
     .           nesum    =   2                     ,
     .           necoil   = 122                     ,
     .           mxtbcmhd = 2000) ! must match keqmax
c
