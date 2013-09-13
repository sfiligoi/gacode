!
!
    MODULE mhdpar
      USE param,only : kj,nw  => nwmhd ,nh =>nhmhd
      PUBLIC
! --- old INCLUDE file mhdpar.i 
!
! --- NOTE: PARAMETER statements also appear in FREYA and NEWBOXX
! --- NOTE: kpsi must be greater than or equal to max of nw and kj
! --- NOTE: kside1,nstark1: for fixed boundary cases mxtbcmhd can
!           be set to a large number ( to accomodate a large list
!           of times in timeqbcd) Consequently the experimental
!           data arrays associated with the free boundary calculations
!           become very large. However  if mxtbcmhd is large enough
!           then we start using virtual memory which destroys the
!           performance of the code to an unacceptable degree.
!           If free boundary equilibria are to be done then
!           set keqmax=mxtbcmhd <~ 300 and set nstark = nstark1
!          set nflxeqbcd = kside, and set magpri  = m
!
!
!  
      integer*4, parameter ::                    &
             nwh      = nw * nh,                 &
             kpsi     =  MAX(nw,kj),             &  ! MAX (nw, kj)
             nfcoil   =  18,                     &
             nrstep   =   1,                     &
             kcoilq   = nfcoil * nfcoil,         &
             kside    =   2 * (nh + nw - 2),     &
             nflxeqbcd = kside,                  &
             nstark1  =   8                     ,&
             nstark   =   nstark1,               & ! use for fixed boundary cases
                nsilop   =  41,                  &
                nvessel  =  24,                  &
                nrogow   =   1,                  &
                magpr322 =  60,                  &
                magpri67 =  29,                  &
                magpri1  = magpr322 + magpri67,  &
                magpri   = magpri1,              &
                magpr2   = magpr322,             &    
                maxobsrv = nsilop + magpr2 + nfcoil, &
                nesum    =   2,                      &
                necoil   = 122,                      &
                mxtbcmhd = 2000                   ! must match keqmax
!
    END MODULE mhdpar
