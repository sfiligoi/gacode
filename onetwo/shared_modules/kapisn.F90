       MODULE kapisn_module
  
         USE nrtype,                            ONLY: DP,I4B
         
         CONTAINS 



      SUBROUTINE kapisn_driver(j,jneo)
!------------------------------------------------------------------------------
! Original from jek 17-sept-99 version 1.0
!
!... Calculates neoclassical chii from model predicted profiles
!... This adapted version of subroutine KAPISN calculates chii at
!... one position instead of the full profile.
! ---  f90 version for Onetwo,gcnmp  HSJ
!    Models:
!      nkimod = 1 Rutherford neoclassical ion conductivity - Julich
!      nkimod = 2 Modified Hazeltine-Hinton as used in BALDUR code
!      nkimod = 3 Bolton model implemented by McCune 7/81
!      nkimod = 4 Chang-Hinton model, Phys Fluids 25, 1982
!      nkimod = 5 Modified Chang-Hinton model for Zeff > 1, ppf 1986
! INPUT from argument list
!   jneo       nkimod is set to jneo -2, jneo in [3,7] is required, not checked for
!   j          zone center grid point index, valid input from j=1 to nj-1
!
! OUTPUT to modules
!    note that j is zone center index:
!    For Onetwo:
!          d(nion+1,nion+1,j)
!          d(nion+2,nion+2,j)
!          xkeneo(j)   
!          xkineo(j)   
!          chiineo(j)       
        
!    For gcnmp
!          dcoefp(nion+1,nion+1,j)
!          dcoefp(nion+2,nion+2,j)
!          iwangrot =  -3  
!---------------------------------------------------------------------------------
!

#ifndef GCNMP  ! for onetwo:

       USE numbrs,                                     ONLY :  nprim,nion,nimp
       USE tcoef,                                      ONLY :  d,xkeneo,xkineo,            &
                                                               chiineo   
       USE tfact,                                      ONLY :  wneo
       USE soln,                                       ONLY :  ene,en,te,ti
       USE ions,                                       ONLY :  atw,atomno,zeff
       USE extra,                                      ONLY :  q
       USE machin,                                     ONLY :  rmajor,btor
       USE ifs,                                        ONLY :  rhod_ifs
       USE mesh,                                       ONLY :  grho2_mesh
       USE neo2d,                                      ONLY :  elong_r
       USE tordlrot,                                   ONLY :  iwangrot
#else          ! for gcnmp:

       USE ions_gcnmp,                                 ONLY :  nprim,nion,nimp,atw,atomno, &   
                                                               zeff
       USE dep_var,                                    ONLY :  en,ene,te,ti
       USE plasma_properties,                          ONLY :  mhd_dat,dischg,dcoefp                    
       USE neo_tport,                                  ONLY :  iwangrot
       USE Vector_class,                               ONLY :  length_Vector
       USE neo_tport,                                  ONLY :  wneo
#endif

               ! common modules :
       USE common_constants,                           ONLY : zeroc





      IMPLICIT NONE

      INTEGER(I4B) jneo,j
      CHARACTER line*132
!
      INTEGER lprint, nkimod_nc, ng_nc, istringer_nc
      REAL(DP) aimp_nc, xzimp_nc, btf_nc,              &
               rminor_nc, rmajor_nc, dshafdr_nc
!
      REAL(DP) rhoel_nco, rhi_nco, rhoi_nco,           &
               te_nco, ti_nco, zeff_nco, q_nco,        &
               xkapi_nco, xnstari_nco, xnstare_nco
!
      REAL(DP)  rhob_exp
!
      REAL(DP) elong_exp, gradrhosq_exp &
       , xkapi_m, xnstari_m &
       , xnstare_m, zfluxlim_nco, chiineol, chieneo



      REAL(DP),DIMENSION(:) :: aplasm_nc(nprim),rhob_nco(nprim),ena(nprim)
      REAL(DP) ensum,enea,zimp,dimp,aimp
      INTEGER(I4B) i,k,jmax

!--------------------------------------------------------
! -- map variables into kapisn names for easier updates
!--------------------------------------------------------
      ng_nc         = nprim
      nkimod_nc     = jneo-2                           ! call with 3 ,=jneo <= 7 only



!     individual ion input:
      ensum = zeroc
      DO i=1,nprim
         ena(i)  = 0.5_DP*(en(j+1,i)+en(j,i))           ! #/cm^3 onetwo, #/m^3 gcnmp
         ensum   = ensum + ena(i)
      ENDDO

      aplasm_nc(1:nprim) = atw(1:nprim)

!     average impuriy:
      zimp=zeroc     ;    dimp = zeroc  ; aimp = zeroc
      DO i=1,nimp
         k =nprim+i
         zimp    = zimp + atomno(k)
         dimp    = dimp + 0.5_DP*(en(j+1,k)+en(j,k))     ! #/cm^3 onetwo, #/m^3 gcnmp
         aimp    = aimp + atw(k)
      ENDDO
      xzimp_nc   = zimp/nimp
      aimp_nc    = aimp/nimp

              




      te_nco     = 500._DP*(ti(j+1)+te(j))        ! ev , zone ctr
      ti_nco     = 500._DP*(ti(j+1)+te(j))        ! ev , zone ctr
      zeff_nco   = (zeff(j+1)+zeff(j))*0.5_DP

 
      dshafdr_nc = zeroc
        
! 
#ifndef GCNMP 
      rmajor_nc= rmajor*100._DP                   ! cm
      IF(rmajor_nc .GT. 700)THEN                  ! ASSUME RMAJOR,RMINOR IN SAME UNITS
         rmajor_nc = rmajor_nc*0.01_DP            ! INPUT WAS CM NOT M
      endif
      rminor_nc         = (rhod_ifs(j+1)+rhod_ifs(j))*50._DP                             !cm
      gradrhosq_exp     = (grho2_mesh(j+1)+grho2_mesh(j))*0.5_DP
      elong_exp         = (elong_r(j+1)+elong_r(j))*0.5_DP
      q_nco             = ABS((q(j+1)+q(j)))*0.5_DP
      btf_nc            = ABS(btor)        
      enea              = (ene(j+1)+ene(j))*0.5_DP                                       ! #/cm^3   
      IF(btf_nc .gt. 15._DP)btf_nc = btf_nc*1.e-4 ! input probably in gauss
#else
       jmax             = length_Vector(dischg%rmajavnj)
       rmajor_nc        = dischg%rmajavnj%data(jmax)*100._DP                             !cm
       rminor_nc        = (dischg%rminavnj%data(j+1)+ dischg%rminavnj%data(j))*50._DP    !cm
       elong_exp        = (dischg%elongxnj%data(j+1)+ dischg%elongxnj%data(j))*0.5_DP
       gradrhosq_exp    = (dischg%grho2nj%data(j+1)+ dischg%grho2nj%data(j))*0.5_DP
       btf_nc           = mhd_dat%btor                                                   !tesla
       q_nco            = ABS((mhd_dat%q_value%data(j+1)+mhd_dat%q_value%data(j)))*0.5_DP
       enea             = (ene(j+1)+ene(j))*0.5d-6                                       ! 1/cm^3
       dimp             = dimp*1.d-6                                                     ! 1/cm^3
       ena(:)           = ena(:)*1.d-6                                                   ! 1/cm^3
       ensum            = ensum*1.d-6                                                    ! 1/cm^3
#endif

      rhi_nco           = dimp/nimp
      rhoel_nco         = enea
      rhoi_nco          = dimp + ensum
      rhob_exp          = rhoi_nco
      rhob_nco(:)       = ena(:)
      istringer_nc      = 0
!




!     write( 126,FMT='("j,rmajor,rminor=",i5,2(x,1pe12.3))')j,rmajor_nc,rminor_nc
!     write( 126,FMT='("btor,gradrhosq_exp,elong =",3(x,1pe12.3))')btf_nc,gradrhosq_exp,elong_exp
 
!
!..initialize output variables
!
      xkapi_m          = 0.D0
      xnstari_m        = 0.D0
      xnstare_m        = 0.D0
      zfluxlim_nco     = 0.D0
      chiineol         = 0.D0
!
!


!
          CALL kapisn(       &
                nkimod_nc,   &      !kapi model desired; default=2 
                aimp_nc,     &      !atomic mass of av. impurity; default 
                xzimp_nc,    &      !atomic number of av. impurity; defau 
                aplasm_nc,   &      !atomic mass of hyd. species; default 
                ng_nc,       &      !number of hyd. species; default=1 
                rhoel_nco,   &      !electron density (cm**-3) 
                rhob_nco,    &      !hyd. spec. den s (cm**-3) 
                rhi_nco,     &      !avg. impurity s density (cm**-3) 
                rhoi_nco,    &      !total ion density (cm**-3) 
                te_nco,      &      !Te (ev) 
                ti_nco,      &      !Ti (ev) 
                zeff_nco,    &      !plasma zeff 
                q_nco,       &      !safety factor 
                btf_nc,      &      !tor field at cntr of outer flux (tes 
                dshafdr_nc,  &      !der. of shaf shift of outer zone bou 
                rminor_nc,   &      !plasma minor radius (cm) 
                rmajor_nc,   &      !major radius (cntr of outer flux) (c 
                istringer_nc,&      !Stringer correction; default=0 
                                    ! Output: 
                xkapi_nco,   &      !o Neoclass ion thermal diffusivity ( 
                xnstari_nco, &      !o nu-star-ions, 
                xnstare_nco, &      !o nu-star-elecs, 
                zfluxlim_nco)       !o flux lim flow max temp grad length
!
!... calculate conductivity, diffusivity, collisionalities
!
         xkapi_m=xkapi_nco*1.D-4*rhob_exp*1.D6/elong_exp      ! 1/(m*s)
         chiineol=xkapi_m/rhob_exp/1.D6/gradrhosq_exp         ! m**2/s
         IF (q_nco .LT. 1) THEN
           chieneo=chiineol                                   ! m**2/s
         ELSE
           chieneo=1.D0/SQRT(aplasm_nc(1))/42.D0*chiineol     ! m**2/s
         ENDIF
!
         xnstare_m=xnstare_nco
         xnstari_m=xnstari_nco


!---------------------------------------------------------------------
! -- map back to onetwo or gcnmp  vars:
!---------------------------------------------------------------------

#ifndef GCNMP
         d(nion+1,nion+1,j)          = wneo(2,2)*1.D4*chieneo*enea      ! electron energy,1/(cm*s)
         d(nion+2,nion+2,j)          = wneo(3,3)*xkapi_m*0.01_DP        ! ion energy,1/(cm*s)
         xkeneo(j)                   = d(nion+1,nion+1,j)
         xkineo(j)                   = d(nion+2,nion+2,j)
         chiineo(j)                  = chiineol*1.D4
#else
         dcoefp(nion+1,nion+1,j)     = wneo(2,2)*chieneo*enea*1.d6      ! 1/(m*s)
         dcoefp(nion+2,nion+2,j)     = wneo(3,3)*xkapi_m                ! 1/(m*s)
#endif

         iwangrot                    =  -3                              ! forces chiwneo  to chiineo in sub ddiffuse





!
!... printout
!
!        WRITE(126,30)
!        IF (nkimod_nc .EQ. 2) WRITE(126,50)
!        IF (nkimod_nc .EQ. 3) WRITE(126,60)
!        IF (nkimod_nc .EQ. 4) WRITE(126,70)
!        IF (nkimod_nc .EQ. 5) WRITE(126,80)
!        WRITE(126,90)
!        WRITE(126,100)
!        WRITE (126,110) rhoel_nco, rhoi_nco  &
!       ,       te_nco, ti_nco, zeff_nco &
!       ,       q_nco, chiineol, chieneo
!
 30     FORMAT('Results from KAPISN model',/ &
       ,       '*************************')
 50     FORMAT('Using modified Hinton-Hazeltine diffusivity',/)
 60     FORMAT('Using C. Bolton diffusivity',/)
 70     FORMAT('Using original Chang-Hinton diffusivity',/)
 80     FORMAT('Using impurity corrected Chang-Hinton diffusivity',/)
 90     FORMAT(t3,'densities [10^19/m^3]',/ &
       ,      t3,'temperatures [keV]',/ &
       ,      t3,'diffusivity [m^2/s]')
 100    FORMAT(/,t126,'ne',t15,'ni',t24,'Te',t33,'Ti' &
       ,      t41,'zeff',t51,'q',t58,'chiineo',t67,'chieneo')
 110    FORMAT (8(1pe9.3,x))
!
 900  CONTINUE
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      END SUBROUTINE kapisn_driver




        SUBROUTINE kapisn( &
                nkimod,    & !model number desired 
                aimp,      & !atomic mass of avg. impurity 
                xzimp,     & !atomic number of avg. impurity 
                aplasm,    & !array of atomic masses of hyd. species 
                ng,        & !number of hydrogen species 
                rhoel,     & !zone centered electron density (cm**-3) 
                rhob,      & !hyd. spec. den_s (cm**-3) 
                rhi,       & !avg. impurity density 
                rhoi,      & !total ion density 
                te,        & !Te (ev) 
                ti,        & !Ti (ev) 
                zeff,      & !plasma zeff 
                q,         & !array of safety factor 
                btf,       & !tf (tesla) fld at cntr of outer flux 
                dshafdr,   & !derivative of shaf shift of outer zone bdr. 
                rminor,    & !plasma minor radius (cm) 
                rmajor,    & !major radius (cntr of outer flux) (cm) 
                istringer, & !stringer correction 
                xkapi,     & !Neoclassical ion thermal diffusivity (cm**2/ 
                xnstari,   & !nu-star-ions, 1 at  banana-plateau transitio
                xnstare,   & !nu-star-elecs,1 at  banana-plateau transitio 
                zfluxlim)    !flux limited max temp gradient length (cm)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  v2.67    CALCULATES NEO-CLASSICAL ION HEAT CONDUCTIVITY
!           TO ACCOUNT FOR SEVERAL HYDROGENIC SPECIES, THE
!           ATOMIC MASS HAS BEEN REDEFINED FOLLOWING
!           RUTHERFORD'S APPROACH USED IN THE BALDUR CODE
!
!... This is the one-gridpoint version of kapisn.
!... Adapted by J. Kinsey July 20, 1999
!
!... Declare variables
!
        IMPLICIT NONE
!
        INTEGER(I4B)  ng, nkimod, istringer, ih
        REAL(DP) aplasm(ng), rhob(ng), &
          rhoel, rhoi, te, ti, q, &
          dshafdr, xkapi, xnstari, &
          xnstare, rhi, zeff, aimp, xzimp, &
          zfluxlim, clogi, btf, rminor, rmajor
!
        REAL(DP) zk20, za2, zb2, zc2, &
          zlamd, aratio, aratii, bpol, zmass, znhy, &
          ztauh, zrhoh, zeps, zhstar, zalpha, zcoef, &
           zkhbp, zkhpc, zkips, zchrg, zchrg4, ztau, &
          zvth, zeps32, znuhat, znstar, zrho, zrho2, &
          rminloc, zbanana, zstring, zkb0b2, zk2st, &
          zfps, zepsc, zb2b0i, zf, zhp, zepsq, zk2ht, &
          zmstar, zk2, znhatb, zkcoef
!
!... Revision History
!
!  REVISED:  JUN-20-79  D. MC CUNE   NEOCLASSICAL CONDUCTIVITY
!          IS MULTIPLIED BY NAMELIST-INPUT FUDGE FACTOR XKFAC
!          AFTER IT IS CALCULATED. (XKFAC DEFAULTS TO 1.0)
!
!  REVISED:  JAN-6-81  D. MC CUNE  *** NEW COMMON VARIABLE ***
!         "NKIMOD" TO SELECT ION KAPA MODEL:  NKIMOD=1==>
!         OLD MODEL (RUTHERFORD--JULICH); NKIMOD=2==> MODIFIED
!         HAZELTINE-HINTON, AS USED IN BALDUR CODE.  THE
!         HAZELTINE-HINTON FORMULA FIT TO NEOCLASSICAL ION
!         THERMAL CONDUCTIVITY INTEGRALS WAS INTENDED FOR A
!         SINGLE SPECIES PLASMA; HERE IT IS USED APPROXIMATELY
!         GENERALIZED FOR MULTI-SPECIES PLASMA, FOLLOWING THE
!         APPROACH IN BALDUR; VALID ONLY IF ZEFF IS CLOSE TO
!         1.0.  (NOTE THESE COMMENTS ALSO APPLY TO JULICH FIT;
!         FOR SINGLE SPECIES ZEFF=1 PLASMA THE HAZELTINE-HINTON
!         FIT IS SUPERIOR).  THE RANGE OF VALIDITY OF THE
!         APPROXIMATIONS USED HEREIN IS A TOPIC WORTHY OF
!         INVESTIGATION.
!            NOTE:  DEUTERIUM-HYDROGEN AND TRITIUM-DEUTERIUM
!         MIXING MAY ALSO EFFECT VALIDITY OF HAZELTINE-HINTON
!         (AND JULICH) APPROXIMATIONS
!  10-JULY 1981  D. MC CUNE
!    ADDING NKIMOD=3 OPTION TO USE C.W. BOLTON NEOCLASSICAL
!  CONDUCTIVITY MODEL: SEE MAY 1981 U. OF TEXAS DISSERTATION
!  "THE ION THERMAL CONDUCTIVITY FOR A PURE TOKAMAK PLASMA"
!  BY C.W. BOLTON FOR HIS PHD.
!  ACCESSED BY SETTING INPUT VARIABLE NKIMOD=3.
!  NOTE BALDUR STYLE ZEFF APPROXIMATE- MODIFICATION IN COLLISIONALITY
!  IS STILL APPLIED AS AN ATTEMPT TO DEAL WITH IMPURE PLASMA, BUT
!  THIS IS APPARENTLY INVALID IN THE CONTEXT OF THE BOLTON
!  CALCULATION, THEREFORE NKIMOD=3 SHOULD BE USED ONLY IF ZEFF IS
!  SMALL *******
!
!  REV 10 SEP 1981  D. MC CUNE:
!   NEW PARAMETER XKFAC2 ADDED.  THIS ANOMOLOUS MULTIPLIER IS APPLIED
!  TO KAPA(I) AT EACH SURFACE WHERE Q .LE. 1.
!
!   stan kaye -- adding chang-hinton kapa(i) (NKIMOD=4)
!     hazeltine hinton with shaf shift and finite epsilon corrections
!   PHYS. FLUIDS VOL. 25 SEPT 1982
!
!   DMC - NKIMOD=5:  CHANG HINTON WITH ADDED CORRECTIONS FOR ZEFF.NE.1
!    AS PER LETTER 20 JULY 1984 FROM DR CHANG TO DR GOLDSTON
!
!       KI*NI = (KH(BP) + KH(PS + CL)) * NH + KI(PS) * NIMP
!         (not sure what above comment refers to - RJG 9/7/87)
!
!   NKIMOD=5 model modified to agree with Chang and Hinton Physics
!       of Fluids article, (1986) p 3314
!       RJG 9/7/87
!
!-----------------------------------------------------------------
!       Modifications to make this a viable subroutine for SNAP
!       Eliminated anomalous factors in this code (all external now)
!       Note TRANSP version defined aratio as r/R, this one is
!       corrected in notation, but gives same physics
!       Rob Goldston 8/25/84
!
!       Modified to include Zeff(r).
!       Mike Zarnstorff   4/29/85
!-----------------------------------------------------------------
!
!  HAZELTINE-HINTON COEFFICIENTS
!
        DATA ZK20/.66D0/
        DATA ZA2/1.03D0/
        DATA ZB2/0.31D0/
        DATA ZC2/0.74D0/
!
!-----------------------------------------------------------------
!
!       CALCULATE THE ION COULOMB LOGARITHM
!
        ZLAMD=1.54D10*TI/ZEFF*DSQRT(TE/RHOEL)
!       BOMB SHELTER: RJG 3/26/85
        IF(ZLAMD .GT. 0.0) THEN
            CLOGI=DLOG(ZLAMD)
        ELSE
            CLOGI=17.D0
        ENDIF
!
!       CALCULATE THE NEOCLASSICAL ION-HEAT CONDUCTIVITY
!
        aratio = rmajor/MAX(1.D-10,rminor)
        aratii = 1.D0/aratio
        bpol = btf/(aratio*q)
!
!       AVERAGE HYDROGENIC MASS AND DENSITY
!
        ZMASS=0.D0
        ZNHY=0.D0
        DO 210 IH=1,NG
          ZMASS=ZMASS+RHOB(IH)*APLASM(IH)
          ZNHY=ZNHY+RHOB(IH)
 210    CONTINUE
!
        ZMASS=ZMASS/ZNHY
!
!  BRANCH ACCORDING TO CHOICE OF CONDUCTIVITY MODEL
        GO TO (220,240,240,240,240) NKIMOD
!
 220    CONTINUE
!
!   JULICH  RUTHERFORD
!
        ZTAUH=1.473D7*DSQRT(ZMASS)*TI**1.5D0/ZNHY/CLOGI
        ZRHOH=1.445D-2*DSQRT(ZMASS*TI)/Btf
        ZEPS=ARATIO**1.5D0        !suspicious but as in TRANSP
        ZHSTAR=4.897D-14*ZNHY*CLOGI*ARATII* &
             Btf*rmajor*ZEPS/(BPOL*TI*TI)
        ZALPHA=XZIMP*XZIMP*RHI/ZNHY
        ZCOEF=ZRHOH*ZRHOH*Q*Q/ZTAUH
        ZKHBP=ZCOEF*ZEPS*.48D0*(1.D0+1.5D0*ZALPHA)/ &
             (1.D0+.36D0*(1.D0+1.5D0*ZALPHA)*ZHSTAR)
        ZKHPC=ZCOEF*(.7D0/Q/Q+1.13D0+ZALPHA*(.5D0+.55D0/ &
             (.59D0+ZALPHA)))
        ZKIPS=ZCOEF*DSQRT(AIMP/ZMASS)*RHI/ZNHY*ZALPHA*1.13D0
!
        XKAPI=((ZKHBP+ZKHPC)*ZNHY+ZKIPS*RHI)/RHOI
!
!       Calculate nustari as in Hazeltine Hinton anyway

        ZCHRG=1.D0
        ZCHRG4=ZCHRG*ZCHRG*ZCHRG*ZCHRG
        ZTAU=2.083D7*DSQRT(ZMASS*TI*TI*TI)/(ZCHRG4*CLOGI)
!  IMPURITIES CORRECTION  (NORMAL HAZELTINE-HINTON WOULD HAVE
!     HYDROGENIC DENSITY IN DENOMINATOR)
        ZTAU=ZTAU/(ZEFF*RHOEL)
        ZVTH=DSQRT(3.204D-12*TI/(1.67D-24*ZMASS))
        ZEPS=ARATII
        ZEPS32=DSQRT(ZEPS*ZEPS*ZEPS)
!
!  NU-STAR
        ZNUHAT=1.414214D0*ARATII*BTF*RMAJOR/(BPOL*ZVTH*ZTAU)
        ZNSTAR=ZNUHAT/ZEPS32
        xnstari=znstar
        xnstare=znstar*1.414214D0*(TI/TE/zchrg)**2
!
        GO TO 200
!
!  HAZELTINE HINTON / BOLTON / chang hinton
!
 240    CONTINUE
!
        ZCHRG=1.D0
        ZCHRG4=ZCHRG*ZCHRG*ZCHRG*ZCHRG
        ZTAU=2.083D7*DSQRT(ZMASS*TI*TI*TI)/(ZCHRG4*CLOGI)
!  IMPURITIES CORRECTION  (NORMAL HAZELTINE-HINTON WOULD HAVE
!     HYDROGENIC DENSITY IN DENOMINATOR)
        ZTAU=ZTAU/(ZEFF*RHOEL)
        ZVTH=DSQRT(3.204D-12*TI/(1.67D-24*ZMASS))
        ZEPS=ARATII
        ZEPS32=DSQRT(ZEPS*ZEPS*ZEPS)
!
!  NU-STAR
        ZNUHAT=1.414214D0*ARATII*BTF*RMAJOR/(BPOL*ZVTH*ZTAU)
        ZNSTAR=ZNUHAT/ZEPS32
        xnstari=znstar
        xnstare=znstar*1.414214D0*(TI/TE/zchrg)**2
!
        ZRHO=1.445D-2*DSQRT(ZMASS*TI)/(ZCHRG*BPOL)
        ZRHO2=ZRHO*ZRHO

! begin ********************************************
! 8.9.95 fix for  Stringer and flux limited flow
!
! istringer=0 recovers old formula
! istringer=1 accounts for reduction of transport inside
! a banana width of r=0
! zfluxlim is an output giving maximum allowed
! temperature scale length -d ln Ti /d rho < 1/zfluxlim

        rminloc=rmajor*aratii+1.D-10
        zbanana=dsqrt(rminloc/rmajor)*ABS(zrho)
        zstring=(zbanana**.66666D0)*(rminloc**.333D0)
        IF (istringer.EQ.0) THEN
         zeps=rminloc/rmajor
        ENDIF
        IF (istringer.EQ.1) THEN
         zeps=rminloc/rmajor
         IF(rminloc.LE.zbanana)  zeps=zstring/rmajor
        ENDIF
!
! now repeat calc for
!
        ZEPS32=DSQRT(ZEPS*ZEPS*ZEPS)
!
!  NU-STAR
        ZNUHAT=1.414214D0*ARATII*BTF*RMAJOR/(BPOL*ZVTH*ZTAU)
        ZNSTAR=ZNUHAT/ZEPS32
        xnstari=znstar
        xnstare=znstar*1.414214D0*(Ti/Te/zchrg)**2
!       write(126,*) 'xnstari=',xnstari
!       write(126,*) 'bt=',BTF
!       write(126,*) 'bpol=',BPOL
!       write(126,*) 'Ti=',TI
!       write(126,*) 'zmass=',ZMASS
!       write(126,*) 'zvth=',ZVTH
!       write(126,*) 'zeff=',ZEFF
!       write(126,*) 'rhoel=',RHOEL
!       write(126,*) 'ztau=',ZTAU
!       write(126,*) 'r/R=',ARATII
!
! end ***********************************************
!
!  HAZELTINE-HINTON / chang-hinton:
        IF(NKIMOD.EQ.4) THEN
!  CHANG HINTON ORIGINAL VSN:
!
         zkb0b2=(1.D0+1.5D0*(zeps**2+zeps*dshafdr)+ &
                3.D0*zeps**3*dshafdr/8.D0)/ &
               (1.D0+0.5D0*zeps*dshafdr)
         zk2st=(1.D0+2.85D0*dsqrt(zeps)-2.33D0*zeps)*zkb0b2
          ZFPS=ZEPS32
         ENDIF
        IF(NKIMOD.EQ.5) THEN
!  CHANG HINTON WITH IMPROVED IMPURITY CORRECTION AS PER LETTER TO
!  R GOLDSTON 20 JULY 1984 FROM DR C.S. CHANG
!  THE ONLY IMPURITY CORRECTION IN THE ORIGINAL VSN WAS THE ENHANCEMENT
!  TO COLLISIONALITY AND ZTAU=THE ION-ION COLLISION TIME
!
!  Modified to follow Physics of Fluids (1986) p 3314 article by Chang
!       and Hinton, rather than private communication.
!
         zkb0b2=(1.D0+1.5D0*(zeps**2+zeps*dshafdr)+ &
                3.D0*zeps**3*dshafdr/8.D0)/ &
               (1.D0+0.5D0*zeps*dshafdr)
         ZEPSC=DSQRT(1.D0-ZEPS*ZEPS)
         ZB2B0I=ZEPSC*(1.D0+0.5D0*ZEPS*dshafdr)/ &
               (1.D0+dshafdr/ZEPS*(ZEPSC-1.D0))
          ZF=0.5D0*(ZKB0B2-ZB2B0I)/DSQRT(ZEPS)
!         ZEFFAC=ZEFF*RHOEL/RHOI
!         Looks like letter had a different impurity parameter
          zalpha = xzimp**2*rhi/znhy      !as for Rutherford/julich
!  PFIRSCH-SCHLUTER CORRECTION ("HP*F" IN THE LETTER):
!         ZHP=1.D0+(1.D0-1.D0/zalpha)*(1.03D0/(zalpha-0.87D0)-0.68D0)
!         New version of Hp:
          zhp = 1.D0 + 1.33D0*zalpha*(1.D0+0.6D0*zalpha)/ &
                (1.D0+1.79D0*zalpha)
          ZFPS=ZHP*ZF
!  BANANA CORRECTION ("HB/0.66" IN THE LETTER)-- AGREES WITH ORIGINAL
!  VSN IN LIMIT ZEFFAC=1.0
!        zk2st=(0.66D0+(1.88D0*dsqrt(zeps)-1.54D0*zeps)*
!              (0.59D0+0.41D0/zalpha))
!     >         *zkb0b2/0.66D0
!        The guy above needs upgrading:
         zepsq = dsqrt(aratii)
         zepsq= dsqrt(zeps)
         zk2ht = ( 0.66D0*(1.D0+1.54D0*zalpha) + &
                 (1.88D0*zepsq-1.54D0*zeps) &
                 *(1.D0+3.75D0*zalpha) ) * zkb0b2
!        Chang-Hinton 1986 introduces a mu*star, and leaves tau as
!        the bulk-bulk collision frequency
         ztau = ztau * (zeff*rhoel) / znhy
         zmstar = znstar * (1.D0 + 1.54D0*zalpha)*znhy / (zeff*rhoel)
         ENDIF
!
!  ORIGINAL HAZELTINE HINTON-- NO ADDL CORRECTIVE FACTORS:
        IF (nkimod .EQ. 2) THEN
          zk2st=1.D0
          ZFPS=ZEPS32
        ENDIF
!
!  HH AND DESCENDENT MODELS PUT TOGETHER:
        IF((NKIMOD.EQ.2).OR.(NKIMOD.EQ.4)) THEN
          ZK2=ZK20*(zk2st/(1.D0+ZA2*DSQRT(ZNSTAR)+ZB2*ZNSTAR) + ZFPS* &
          ZEPS32*(ZC2*ZC2/ZB2)*ZNSTAR/(1.D0+ZC2*ZEPS32*ZNSTAR) )
          GO TO 245
        ENDIF
        IF(NKIMOD.EQ.5) THEN
          ZK2=ZK20*(zk2ht/0.66D0/(1.D0+ZA2*DSQRT(ZmSTAR)+ZB2*ZmSTAR)+ &
          ZFPS*ZEPS32*(ZC2*ZC2/ZB2)*ZmSTAR/(1.D0+ZC2*ZEPS32*ZmSTAR) )
          GO TO 245
        ENDIF
!
!  BOLTON
        ZNHATB=.75D0*DSQRT(3.1415926D0)*ZNUHAT
!  BOLTON_S COLLISIONALITY DIFFERS FROM HAZELTINE-HINTON_S BY FACTOR
!  3*SQRT(PI)/4 ...
        ZK2=ZKBOLT(ZNHATB,ZEPS)
!
 245    ZKCOEF=ZRHO2/ZTAU
        XKAPI=ZK2*DSQRT(ARATII)*ZKCOEF
!
!*******************************************************
        xkapi=zk2/ztau* &
           dsqrt(zeps)*(zbanana+1.D-10)**2/(zeps+1.D-10)
        IF (istringer.EQ.1..AND.rminloc.LE.zbanana) THEN
         xkapi=xkapi*rminloc/(zbanana+1.D-10)
        ENDIF
        zfluxlim=zbanana+1.D-10
        IF (rminloc.LE.zbanana) zfluxlim=zstring+1.0D-10
!*******************************************************
!
        GO TO 200
!
!
200     CONTINUE
!
        RETURN
!
!
        END SUBROUTINE kapisn
!
!*******************************************************
!
!  BOLTONS KI(NU-HAT,EPS) FUNCTION (THESIS PP 74--78)
!
        REAL(DP) FUNCTION zkbolt(znuhat,zeps)
!
        IMPLICIT NONE
!
        REAL(DP) ZNUHAT, ZEPS, ZSE, ZEPS32, ZTANC, &
          ZTBNC, ZTCNC, ZTDNC, ZNH, ZKINC, &
          ZTAPS, ZTBPS, ZTCPS, ZTDPS, ZKIPS

!  N-C COEFFICIENTS FROM TABLE VI
        REAL(DP) ZANC(3),ZBNC(3),ZCNC(3),ZDNC(3)
!  P-S COEFFICIENTS FROM TABLE VI
        REAL(DP) ZAPS(3),ZBPS(3),ZCPS(3),ZDPS(3)
!
        DATA ZANC/2.441D0,-3.87D0,2.19D0/
        DATA ZBNC/.9362D0,-3.109D0,4.087D0/
        DATA ZCNC/.241D0,3.40D0,-2.54D0/
        DATA ZDNC/.2664D0,-.352D0,.44D0/
!
        DATA ZAPS/.364D0,-2.76D0,2.21D0/
        DATA ZBPS/.553D0,2.41D0,-3.42D0/
        DATA ZCPS/1.18D0,.292D0,1.07D0/
        DATA ZDPS/.0188D0,.180D0,-.127D0/
!
        ZSE=DSQRT(ZEPS)
        ZEPS32=ZSE*ZSE*ZSE
        ZTANC=.66D0+((ZANC(3)*ZSE+ZANC(2))*ZSE+ZANC(1))*ZSE
        ZTBNC=((ZBNC(3)*ZEPS+ZBNC(2))*ZEPS+ZBNC(1))/(ZEPS**.75D0)
        ZTCNC=((ZCNC(3)*ZEPS+ZCNC(2))*ZEPS+ZCNC(1))/ZEPS32
        ZTDNC=((ZDNC(3)*ZEPS+ZDNC(2))*ZEPS+ZDNC(1))/ZEPS32
!
!  N-C PART
!
        ZNH=DSQRT(ZNUHAT)
        ZKINC=ZTANC/ &
        (1.D0+ZTBNC*ZNH+ZTCNC*ZNUHAT+ZTDNC*ZNUHAT*ZNUHAT)
!
        ZTAPS=((ZAPS(3)*ZSE+ZAPS(2))*ZSE+ZAPS(1))*ZEPS32
        ZTBPS=((ZBPS(3)*ZSE+ZBPS(2))*ZSE+ZBPS(1))*ZEPS32
        ZTCPS=((ZCPS(3)*ZSE+ZCPS(2))*ZSE+ZCPS(1))
        ZTDPS=((ZDPS(3)*ZSE+ZDPS(2))*ZSE+ZDPS(1))
!
!  P-S PART
!
        ZKIPS=1.57D0*ZEPS32 + (ZTAPS+ZTBPS*ZNH)/ &
                (1.D0+ZTCPS*ZNUHAT**1.5D0+ZTDPS*ZNUHAT**2.5D0)
!
        ZKBOLT=ZKINC+ZKIPS
!
      RETURN
      END FUNCTION zkbolt


      END MODULE kapisn_module
