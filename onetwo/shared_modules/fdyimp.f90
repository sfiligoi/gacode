     MODULE fdyimp
      USE nrtype ,               ONLY : DP, I4B
      USE common_constants,      ONLY : zeroc,izero

      USE io_gcnmp,               ONLY : ncrt,nout
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
      REAL(DP)                                                       &
                      adot, bdot, ascrip, bscrip,                    &
                      adotmult, f2d1mult,f2d2mult,f2d3mult,          &
                      e2dmult,    qi2dmult,  dlnhdtmult,             &
                      dlnhdrmult,  dlnrdtmult
      REAL(DP), dimension(:),allocatable ::                          &
                      fday2d1, fday2d2, fday2d3 ,dscrip,afdayi,      &
                      bfdayi,cfdayi,dfdt,dgdt,dhdt,dr2idt,drdt,      &
                      drfghdt

!
      CONTAINS 
        SUBROUTINE  allocate_fdy_arrays
          USE grid_class,                           ONLY : nj
            IF(nj .LE. izero)THEN
               write(ncrt,1)nj
               write(nout,1)nj
1              Format(2x,'Error in allocate_fdy_arrays,nj=',i5)
!               CALL STOP('allocate_fdy_arrays',1)
!              STOP not compatible with GCNMP
            ENDIF
            IF(.NOT. ALLOCATED(fday2d1))                                   &
                ALLOCATE(fday2d1(nj),fday2d2(nj),fday2d3(nj),dscrip(nj))
            IF(.NOT. ALLOCATED(afdayi))                                    &
                ALLOCATE(afdayi(nj),bfdayi(nj),cfdayi(nj))
            IF(.NOT. ALLOCATED(dfdt))                                      &
                ALLOCATE(dfdt(nj),dgdt(nj),dhdt(nj))
            IF(.NOT. ALLOCATED(drdt))ALLOCATE(drdt(nj))
            IF(.NOT. ALLOCATED(dr2idt))ALLOCATE(dr2idt(nj))
            IF(.NOT. ALLOCATED(drfghdt))ALLOCATE(drfghdt(nj))

            fday2d1(:)  = zeroc   ; fday2d2(:) = zeroc
            fday2d3(:)  = zeroc   ; dscrip(:)  = zeroc
            afdayi(:)   = zeroc   ; bfdayi(:)  = zeroc
            cfdayi(:)   = zeroc   ; dfdt(:)    = zeroc
            dgdt(:)     = zeroc   ; dhdt(:)    = zeroc
            drdt(:)     = zeroc   ; dr2idt(:)  = zeroc
            drfghdt(:)  = zeroc


        END SUBROUTINE  allocate_fdy_arrays

      END MODULE fdyimp
