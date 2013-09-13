

  MODULE oloss_dat
	USE PARAM,                            ONLY : kf,kj
! --- common block for subroutine OLOSSA ------------------ HSJ 10/17/94
! --- note that npitch is also defined and used in NUBPLT
!
      parameter (npitch = 90)        ! number of pitch angles in table
      parameter (nrxpt  = 25)        ! size of vectors rzxpt and psizxpt
!                                     (the major radius and psi values
!                                      for the radial vectors in olossa)
!
! --- note that bp0 is loaded only if iborb=2 (=> orbit loss calcs done)
!
      REAL*8            pitchlos (npitch,kf), pitchang (npitch   ),       &
                       pitcheng (npitch   ),                              &
                       pitchrmaj(kf),rzxpt(nrxpt), psizxpt(nrxpt),        &
                       bpeq  (kf), bp0(kj),bpoltp(kf), losscalc
!

  END MODULE oloss_dat
