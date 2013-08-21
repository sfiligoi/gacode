     MODULE tmpcom
     USE param,    ONLY : kj
     IMPLICIT NONE
!
! --- ascrip = 1/(rho at limiter at time t)
! --- adot   = d(ascrip)/dt
! --- bscrip = 0   bot = d(bscrip)/dt = 0.0 in versions of code
! --- since about 1982 (bscrip was used for Doublet-type plasmas)
! --- dscrip is defined as -rho*(adot+bdot*rho)/(ascrip+2*bscrip*rho)
! ---                   if the adaptive grid calcs are not selected.
! ---                   If the adaptive grid calcs are selected then
! ---                   dscrip(j) is set equal to drhodt_adaptive(j)
! ---                   (which is calculated from the addative rho grid)
! --- dscrip is the speed of a NORMALIZED flux surface relative to the
! --- magnetic axis.  In time-independent MHD runs the value of rho
! --- at the plasma boundary (i.e., rholim) does not change with time and
! --- all of these parameters go to zero.
! --- fday2d1 =  (bp0/hcap*r)(d/dt)(ln(fcap*gcap*hcap*r))
! --- fday2d2 = -(bp0/hcap*r)(d/dr)(dscrip)
! --- fday2d3 =    (1/hcap*r)(d/dr)(dscrip*bp0)
! --- adotmult a multiplier for adot.  used to gauge effect of this
! --- term on solution (defaulted to 1.0).
! --- similarly f2d1mult,    f2d2mult, and    f2d3mult are multipliers
! --- for    fday2d1    , fday2d2    , and fday2d3  respectively.
!
      REAL*8                                                       &
                      adot, bdot, ascrip, bscrip, dscrip(kj),      &
                      adotmult,                                    &
                      fday2d1(kj), fday2d2(kj), fday2d3(kj),       &
                      f2d1mult,    f2d2mult,    f2d3mult,          &
                      e2dmult,    qi2dmult,  dlnhdtmult,           &
                      dlnhdrmult,  dlnrdtmult
!
      END MODULE tmpcom
