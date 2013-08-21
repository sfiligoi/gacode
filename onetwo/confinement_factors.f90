
  MODULE cfactrs
    USE nrtype,                         ONLY : DP,I4B
    USE numbrs,                         ONLY : nj,nprim
    USE ions,                           ONLY : atw
    USE soln,                           ONLY : en
    USE machin,                         ONLY : rmajor,rminor,elong,btor,   &
                                               kappa
    USE bd_condtn,                      ONLY : totcur
    USE extra,                          ONLY : enebar,etot,eitot,eetot
    USe yoka,                           ONLY : pradt,pbe,pbi,prf
    USE fusion,                         ONLY : pfusdt,pfusdd,pfustt

    REAL(DP) H89pm,H_ITER98y2,H_Petty

    CONTAINS 

      SUBROUTINE confinement_factrs
! -------------------------------------------------------------
! ---
! ----------------------------------------------09/11/06-----HSJ


        IMPLICIT NONE

        REAL(DP) IP,R0,mavg,hmass,hden,hkappa,hcur,hrad,     &
                 hrmaj,hbtor,hpow,htau89p,nbar,Bt0,aratio,   &
                 wtot,wth,Pbrems,Pfus,Rminr,Paux,Ptransport, &
                 Palpha,taue_net,taue_tot,ITER98Y2,TauPetty
        INTEGER(I4B) j,i

        !density weighted avg mass no of primary ions :
        mavg = 0.0
        do j=1,nj
          hmass   = 0.0
          hden    = 0.0
          do i=1,nprim
            hmass = hmass + atw(i) * en(j,i)
            hden  = hden  +          en(j,i)
          end do
          mavg = mavg + hmass/hden
        end do
        mavg = mavg/nj
        IP          = totcur(1)/1.0e6                    ! MA
        R0          = rmajor/100.                        ! m
        Rminr       = rminor/100.                        ! m
        Bt0         = ABS(btor)/1.e4                     ! Tesla 
        aratio      = R0/Rminr                           ! R/a 
        !count d(t,n),d(d,n),t(t,2n) reactions:
        Pfus        = (pfusdt + pfusdd + pfustt)/1.e6    ! MW
        Paux        = pbe + pbi + prf                    ! MW
                                                         ! NOTE: pbe,pbi,prf  are in MW 
        Pbrems      = pradt/1.e6
        Palpha      = (3.5/17.6)*pfusdt/1.e6
        Ptransport  = Paux + Palpha - Pbrems             ! paux +palpha -Pbrems
        wtot        = etot/1.e6                          ! Mj
        wth         = (eetot + eitot)/1.e6               ! MJ
        taue_net    = wth/Ptransport 
        nbar        = enebar/ 1.0e14                     !(10**20/m**3) electron line avg density
        IF(Ptransport .GT. 0)THEN
           ITER98Y2    = 0.0562*(IP**0.93)*(Bt0**0.15)*(Ptransport**(-0.69)) &
                *((10.*nbar)**0.41)*(mavg**0.19)*(R0**1.97)*         &
                (aratio**(-0.58))*(kappa**0.78)



           TauPetty    = (0.028*IP**0.83)*(Bt0**0.07)*((10.*nbar)**0.49)     &
                *(Ptransport**(-0.55))*(mavg**0.14)*(R0**2.11)*      &
                (aratio**(-0.30))*(kappa**0.75)

           H_ITER98y2  = taue_net/ITER98Y2
           H_Petty     = taue_net/TauPetty

           hmass       = mavg**.5
           hden        = nbar**0.1                             ! 10**20/m**3
           hkappa      = kappa**0.5
           hcur        = IP **0.85                             ! in MA
           hrad        = Rminr**0.3                            ! in meters
           hrmaj       = R0**1.2                               ! in meters
           hbtor       = Bt0**0.2                              ! in Tesla
           hpow        = Ptransport**(-0.5)
           htau89p     = 0.048*hmass*hcur*hrmaj*hrad*hkappa*hden*hbtor*hpow
           taue_tot    = wtot/Ptransport 
           H89pm       = taue_tot/ htau89p
        ELSE ! ptransport is negative, cant define these quantitites
           H89pm       = 0.0
           H_Petty     = 0.0
           H_ITER98y2  = 0.0   
        ENDIF
!        print *,'mavg,nbar,kappa,IP =',mavg,nbar,kappa,IP
!        print *,'Rminr,R0,bt0,ptransport =',Rminr,R0,bt0,ptransport
!        print *,'wth,taue_net,htau89p,aratio =',wth,taue_net,htau89p,aratio
!        print *,'Pbrems,pfus,palpha =',Pbrems,pfus,palpha
!        print *,'wth,Ptransport,enebar =',wth,Ptransport,enebar
!        print *,'etot,eetot,eitot =',etot,eetot,eitot
!        print *,'pfus,paux,palpha,nbar =',pfus,paux,palpha,nbar
      RETURN
      END SUBROUTINE confinement_factrs

  END MODULE cfactrs
