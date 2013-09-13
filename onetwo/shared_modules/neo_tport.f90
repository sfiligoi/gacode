   MODULE neo_tport  

     USE nrtype,                                     ONLY : Dp,I4B
     USE common_constants,                           ONLY : zeroc,izero
     USE solcon_gcnmp,                               ONLY : itran_max

     IMPLICIT NONE
     REAL(DP)  wneo1,wneo2, wneo3, wneo4,wneo5,wneot
     REAL(DP), DIMENSION(itran_max,itran_max) ::  wneo,neocl_mult
     LOGICAL use_forcebal_chie,use_forcebal_chii
!-------------------------------------------------------------------------------------------
! -- The following variables are used in sub diffuse_coeff and Ch_neo_single_grid_point
! -- They are defined here to get a common point of definition for both routines.
      REAL(DP)   enea,tea,tia,curdna,dte,dti,zeffa,angrota,zin,enasum,dentia,drbp,    &
                 drangrot,epsa,delta,rdelta,delta2,xhm2a,xft,xfc,del32,del32i,        &
                 del52,r2capa,grsq,bp1,bp2,rbpa,safac,one_hundredth,                  &
                 tejou,teaerg,vthe,rhoe,xlam,xnue,tijou,tiserg,xnubyv,xmassi,xmassp,  &
                 cesu,charg4,sumdz,etaspitzer,fcapa,gcapa,hcapa,charge,x1,x4, &
                 xidum,puny,bpa,xmasse,xipsa,tiaerg,xkap11
       REAL(DP),ALLOCATABLE,SAVE,DIMENSION(:) :: ftfc,xnusi,vth,dwifs,ftncl,za,zsqa,  &
                                                 rzsqa,dzdtea,ena,vionzgrd,       &
                                                 rho_ion,xnuse, xnu,ftrap
       REAL(DP),ALLOCATABLE,SAVE,DIMENSION(:,:) :: etaim,xnus
       INTEGER(I4B) iwangrot,ifus,jneo
       REAL(DP)   xk0(3,3), xk01(3,3), xk0i(3,3),                                     &
                 a(3,3), a1(3,3), a4(3,3), b(3,3), b1(3,3), b4(3,3),                  &
                 c(3,3), c1(3,3), c4(3,3)
       REAL(DP)   xk(3,3), xia(3,3)
       DATA xk01/1.04,0.  ,0.  ,1.20,2.55,0.  ,2.30,4.19,1.83/, &
           xk0i/0.73,0.  ,0.  ,0.73,1.46,0.  ,1.46,3.65,1.46/, &
           a1  /2.01,0.  ,0.  ,0.76,0.45,0.  ,1.02,0.57,0.68/, &
           a4  /2.30,0.  ,0.  ,0.80,0.46,0.  ,0.79,0.48,0.47/, &
           b1  /1.53,0.  ,0.  ,0.67,0.34,0.  ,0.75,0.38,0.32/, &
           b4  /0.98,0.  ,0.  ,0.42,0.22,0.  ,0.56,0.33,0.20/, &
           c1  /0.89,0.  ,0.  ,0.56,0.43,0.  ,1.07,0.61,0.66/, &
           c4  /0.74,0.  ,0.  ,0.48,0.30,0.  ,0.51,0.28,0.51/
       DATA  puny   / 1.0e-6      /
! -- end common definition for sub diffuse_coeff and Ch_neo_single_grid_point
!-------------------------------------------------------------------------------------------

     CONTAINS 

     REAL(DP) FUNCTION f1i(zin, x1, xidum) 
       REAL(DP) zin,x1,xidum
       f1i =  xidum + zin * (x1-xidum)
     END FUNCTION f1i

     REAL(DP) FUNCTION  f14(zin, x1, x4)
       REAL(DP) zin,x1, x4
       f14 =  1.33333*(1.0-zin)*x4 +0.333333*(4.0*zin-1.0)*x1
     END FUNCTION f14


      SUBROUTINE neocl_init
       USE solcon_gcnmp,                            ONLY : include_neocl
       IMPLICIT NONE
       REAL(DP) wc
        include_neocl = izero ; wc = zeroc
        wneo(:,:) = 1.0_DP 
        wneo(:,:) = wneo(:,:) * neocl_mult(:,:)
        wc = SUM(wneo)
        IF(wc > zeroc) include_neocl = 1_I4B
        wneo1 = wneo(1,1)+wneo(1,2)+wneo(1,3)+wneo(1,4)+wneo(1,5)
        wneo2 = wneo(2,1)+wneo(2,2)+wneo(2,3)+wneo(2,4)+wneo(2,5)
        wneo3 = wneo(3,1)+wneo(3,2)+wneo(3,3)+wneo(3,4)+wneo(3,5)
        wneo4 = wneo(4,1)+wneo(4,2)+wneo(4,3)+wneo(4,4)+wneo(4,5)
        wneo5 = wneo(5,1)+wneo(5,2)+wneo(5,3)+wneo(5,4)+wneo(5,5)
        wneot = 1._DP
        
       RETURN 

      END SUBROUTINE neocl_init




         SUBROUTINE  CH_neo(j,sumdz,xnuse,xnusi,xnus,ftrap,xk,xhm2a,ena,enasum,     &
                      za,zeffa,wneo,fcapa,hcapa,rbpa,wneo1,wneo2,wneo3,wneo4,wneo5, &
                      dwifs,tia,tea,teaerg,enea,del32,delta,rdelta,xipsa,bpa,xnue,  &
                      rhoe,rho,xkap11,xnu,ifus,zsqa,xia,jboot,iwangrot,delta2,      &
                      del32i)
!  
! ----------------------------------------------------------------------
!
!      Chang- Hinton                  neoclassical transport
! ----------------------------------------------------------------------
!
      USE nrtype,                 ONLY : DP,I4B

      USE plasma_properties,      ONLY : mhd_dat,dischg,xdchitot,xketot,   &
                                         xchietot,xchiitot,xkangrot,dcoefp,&
                                         chiwneo

      USE ions_gcnmp,             ONLY : nion,dmassden,nprim,atw

      USE source_terms_gcnmp,     ONLY : dneo,xkeneo,qieneo,xkineo,chiineo

      USE curden_terms,           ONLY : q,eta

      USE common_constants,       ONLY : pi,Permeability

      USE solcon_gcnmp,           ONLY : ntot,dudr,constd_values

      USE grid_class,             ONLY : ra,nj,r

      USE common_constants,       ONLY : Proton_Mass 

      IMPLICIT NONE
      REAL(DP) xnuse(*),ftrap(*),xk(3,*),ena(*),za(*),wneo(itran_max,*),  &
               dwifs(*),xnusi(*),xnus(nion,*),di(nion),xnu(*),            &
               rho(*),zsqa(*)
      REAL(DP) d11,d12,d13,d22,d23,betape,sumdz,tt,de,x,y,                &
               betag,betagx,xnu23,xhm2a,wneo1,wneo2,wneo3,wneo4,wneo5,    &
               zeffa,fcapa, hcapa,rbpa,tia,tea,sumen,enea,zefcor,         &
               chfact,xk2,xk20,a2,b2,c2,del32,delta,rdelta,xipsa,         &
               bpa,xnue,rhoe,xkap11,zetax,teaerg,enasum,tejou,delta2,     &
               xmuib,xmuip,xmui,del32i, xkangavg,xkangrob
      REAL(DP) xia(3,3),tiny,cesu,xmult,u0
      INTEGER(I4B) j,i,k,l,ifus,jboot,iwangrot,kfar
      DATA xk20,a2,b2,c2/0.66,1.03,0.29,0.74/


      u0    = Permeability
      tejou = teaerg*1.e-7     
      tt    = 1.0 + tia/(zeffa*tea)
      de    = rdelta*rhoe**2*xnue                                  ! m**2/sec
      DO 206 i=1,nion
  206 di(i) = rdelta*rho(i)**2*xnu(i)
      x     = SQRT (xnusi(j))
      y     = xnusi(j)**2*delta**3
      betag = ((1.17-0.35*x)/(1.0 + 0.7*x)-2.1*y)/(1.0 + y)
      xnu23 = 1.0 + xnuse(j)**2*delta**3
      betagx = betag/xnu23
      zetax  = 1.0/(xkap11*(1.0-ftrap(j)))
      betape = enea*tejou*2.*u0/bpa**2
!
!  calculate dmn factors
!
      d11 = xk(1,1) + 0.5 * zetax*rdelta*xk(1,3)**2/xhm2a
      d12 = xk(1,2) + 0.5 * zetax*rdelta*xk(1,3)*xk(2,3)/xhm2a
      d22 = xk(2,2) + 0.5 * zetax*rdelta*xk(2,3)**2/xhm2a
      d13 = zetax*xk(1,3)/(xhm2a*betape)
      d23 = zetax*xk(2,3)/(xhm2a*betape)
!
! calculate the transport coefficients for ion particle transport.
!
!     obtain terms arising from electron-ion collisions.
!     USE expressions due to Rawls, Chu, and Hinton.
!
      IF (wneo1 .EQ. 0.0)  go to 240
       DO  i=1,nion
          x = ena(i)*de*za(i)/zeffa
          dcoefp(i,i,j) =       dcoefp(i,i,j) + wneo(1,1)*x*tia*d11/(za(i)*tea*ena(i))
          DO k=1,nion
             xmult = 1.0_DP
             IF (ifus .EQ. -1 .AND. k .NE. i)  xmult = 0.0
             dcoefp(i,k,j) = dcoefp(i,k,j) + xmult * wneo(1,1) * x * za(k) * d11 / enea
          END DO
          !
          k        = nion + 1
          dcoefp(i,k,j) = dcoefp(i,k,j) + wneo(1,2)*x*((-1.5+sumdz)*d11+d12)/tea
          k        = nion + 2
          dcoefp(i,k,j) = dcoefp(i,k,j) + wneo(1,3)*x*(1.0-betagx)*d11/(za(i)*tea)
          k        = nion + 3
          dcoefp(i,k,j) = dcoefp(i,k,j) + wneo(1,4)*x*d13/(fcapa**2*hcapa*rbpa)
      ENDDO
!
!     SAVE diffusion coefficient of first ion species
!
      dneo(j) = dcoefp(1,1,j)
  240 CONTINUE
!
! calculate the transport coefficients for electron energy transport.
!
!     obtain terms arising from electron-ion collisions.
!     USE expressions due to Rawls, Chu, and Hinton, except that
!     convection is subtracted out.
!
      IF (wneo2 .EQ. 0.0)  go to 440
      i = nion+1
      x = enea*tea*de
      DO k=1,nion
          dcoefp(i,k,j) = dcoefp(i,k,j)+ wneo(2,1)*x*tt*za(k)*(d12-2.5*d11)/enea
      ENDDO
      k        = nion+1
      dcoefp(i,k,j) = dcoefp(i,k,j)+ wneo(2,2)*x*((-1.5+tt*sumdz)*(d12-2.5*d11) &
                                        +(d22-2.5*d12))/tea
      k        = nion+2
      dcoefp(i,k,j) = dcoefp(i,k,j)+ wneo(2,3)*x*(1.0-betagx)*(d12-2.5*d11) &
                                                  / (zeffa*tea)
      k        = nion+3
      dcoefp(i,k,j) = dcoefp(i,k,j)+ wneo(2,4)*x*(d23-2.5*d13) &
                                      /(fcapa**2*hcapa*rbpa)
!
!     SAVE electron thermal conductivity
!
      dwifs(j)  = dcoefp(i,i,j)
      xkeneo(j) = dcoefp(i,i,j)
!
! calculate the transport coefficients for ion energy transport.
!
!     obtain terms arising from electron-ion collisions.
!     USE expressions due to Rawls, Chu, and Hinton.
!
  440 i = nion + 2                          ! TI equation
      x = -enea*tia*de*betagx/zeffa         ! kev/(m sec) 
      IF (wneo3 .EQ. 0.0 .AND. iwangrot .NE. -3)  go to 640
      sumen = 0.0
      DO  k=1,nion
        sumen = sumen + ena(k)
        dcoefp(i,k,j) = dcoefp(i,k,j) +x*wneo(3,1)*(za(k)+tia/tea)*d11/enea
      ENDDO
      k        = nion + 1
      dcoefp(i,k,j) = dcoefp(i,k,j) +x*wneo(3,2)*((-1.5+sumdz)*d11+d12)/tea
      k        = nion + 2
      dcoefp(i,k,j) = dcoefp(i,k,j) +x*wneo(3,3)*(1.0-betagx)*sumen*d11/(enea*tea)
      k        = nion + 3
      dcoefp(i,k,j) = dcoefp(i,k,j) +x*wneo(3,4)*d13/(fcapa**2*hcapa*rbpa)
!
!     calculate partial ion heat flux for USE in source
!
      DO  k=1,ntot
         qieneo(j) = qieneo(j) - dcoefp(i,k,j)*dudr(k,j)               !kev/(m**2 sec)
      ENDDO
!
!     obtain terms arising from like ion-ion collisions.
!     USE expression due to Chang and Hinton.
!
      zefcor = zeffa*enea/(zsqa(1)**2*ena(1))
      chfact = 1.0 + 2.85*rdelta - 2.33*delta
      xk2    = xk20*((chfact/(1.0 + a2 * SQRT (xnus(1,j)*zefcor) &
              +b2*xnus(1,j)*zefcor)) &
              *xia(1,1)/(1.95*rdelta) &
              +0.5 * c2**2*xnus(1,j)*zefcor*delta*xipsa &
              /(b2*(1.0 + c2*xnus(1,j)*zefcor*del32)))
      DO  l=1,nion
         xkineo(j) = xkineo(j) + xk2*zeffa*enea*di(l)/zsqa(l)**2
      ENDDO
      chiineo(j)= xkineo(j)/enasum             
      xkineo(j) = xkineo(j)*wneo(3,3)

      dcoefp(i,i,j) = dcoefp(i,i,j) + xkineo(j)

      IF(wneo3 .EQ. 0)THEN ! only want chiineo for toroidal rotation in this case
         DO k =1,ntot
            dcoefp(i,k,j) = 0.0_DP
         ENDDO
      ENDIF
  640 CONTINUE
!
      IF (ABS (0.5*(q(j+1)+q(j))) .LT. 1.0    ) &
          dwifs(j) = xkineo(j) ! set dwifs to ion neoclassical value
 
!
!  add classical ion conduction;
!          helium plasmas may go wrong here (among other places)
!

!
! calculate the transport coefficients corresponding to Ohm's law.
!     use expressions due to Rawls, Chu, and Hinton.
!
      kfar = nion + 3
      x = zetax*de*bpa/ra(j)             !     tesla m/sec

      IF (jboot== 0)THEN                 ! jboot = 0 is Hinton+Hazel. bootst
         DO k=1,nion
            dcoefp(kfar,k,j) = dcoefp(kfar,k,j) + x*wneo(4,1)*tt*za(k)*xk(1,3)/(2.0*enea)
         ENDDO
         k        = nion+1
         dcoefp(kfar,k,j) = dcoefp(kfar,k,j) + x*wneo(4,2)*((-1.5+tt*sumdz)*xk(1,3) &
              +xk(2,3))/(2.0*tea)
         k        = nion+2
         dcoefp(kfar,k,j) = dcoefp(kfar,k,j) + x*wneo(4,3)*(1.0-betagx)*xk(1,3) &
              /(2.0*zeffa*tea)
      ENDIF


      !pick up only diagonal term if jboot > 0
      !for Sauter and Nclass above  contributions are handled  
      !in other routines (either as dcoefp terms or as source terms)
      !see 4.2-16 of GAA, units 1/sec ,bpa dependence cancels due to x

      ! eta in ohm m,u0 in H/m
      dcoefp(kfar,kfar,j) = wneo(4,4)* eta(j)/(u0*hcapa*(fcapa*ra(j))**2)        !  1./sec


!
! ----------------------------------------------------------------------
! --- calculate transport coefficients for toroidal rotation.
! --- assume diagonal model for now with neoclassical expression given by
! --- Hinton and Wong and Wong (Phys. Fluids 28,3082(1985) and 30,818(1988)
! --- drift wave model of mattor-diamond (phys fluids 31,1180(1988)
! --- or use empirical take your pick type input model:
! ----------------------------------------------------------------------
!
      wneo5l : IF(wneo5 .GT. 0)THEN
         i        = nion + 4
         dmassden(j) = 0.0
         !
         ! --- determine momentum diffusivity from neoclassical theory
         !
         xmuib      = 0.1 *  delta2
         xmuip      = 0.6 * (bpa/mhd_dat%btor)**2
         chiwneo(j) = 0.0
         !
         !***  do k=1,nion
         DO k=1,nprim
            !
            !       banana regime
            !
            IF      (xnus(k,j) .LE. 1.0) THEN
               xmui = xmuib
            ELSE IF (1.0 .LT. xnus(k,j) .AND. xnus(k,j) .LT. del32i) THEN
               !
               !         plateau regime: assume linear connection (banana to Pfirsch-Schluter)
               !
               xmui = xmuib + (xnus(k,j)-1.0)*(xmuip-xmuib)/(del32i-1.0)
            ELSE             ! Pfirsch-Schluter regime
               xmui = xmuip
            END IF
            xmui       = xmui * rho(k)**2 * xnu(k)
            chiwneo(j) = chiwneo(j)+xmui*atw(k)*Proton_Mass*ena(k)*dischg%rmajor**2 
            dmassden(j)   = dmassden(j)+atw(k)*Proton_Mass*ena(k)
         END DO
         !
         chiwneo(j)   = chiwneo(j)/((dischg%rmajor**2)*dmassden(j)) !cm**2/sec

         !
         ! --- get the diffusion coefficient and the momentum diffusivity
         ! --- based on neoclassical and/or Mattor-Diamond models
         !
         IF (iwangrot .EQ. 0) THEN
            !
            ! ---     neoclassical momentum diffusivity
            !
            xkangrot(j) = chiwneo(j)*wneo(itran_max,itran_max)
            dcoefp(i,i,j)    = dcoefp(i,i,j)+xkangrot(j)*dmassden(j)*dischg%rmajor**2

         END IF
         !
         ! --- xkangrot is set by empirical input model
         !
         IF (iwangrot .NE. 0 .AND. iwangrot .NE. -2) THEN
            xkangavg    = (xkangrot(j)+xkangrot(j+1))*0.5
            dcoefp(i,i,j)    = dmassden(j)*xkangavg*dischg%rmajor**2+dcoefp(i,i,j)
         END IF

         IF(iwangrot .EQ. -3)THEN
            chiwneo(j)   = chiineo(j)        ! m**2/sec
            xkangrot(j)  = chiwneo(j)*wneo(itran_max,itran_max)
            dcoefp(i,i,j) = dmassden(j)*xkangrot(j)*dischg%rmajor**2 ! (kg m)/sec
         ENDIF
         !
         ! --- explicitly zero off-diagonal elements
         !
         DO k=1,i-1
            dcoefp(i,k,j) = 0.0
            dcoefp(k,i,j) = 0.0
         END DO
         !
         ! --- extrapolate to boundary (j = nj, used for printout only)
         !
         IF (j .EQ. nj-1) THEN
            CALL extrap (ra(nj-2), ra(nj-1), r(nj), &
                 xkangrot(nj-2), xkangrot(nj-1), xkangrob)
            xkangrot(nj) = xkangrob
         END IF
     ENDIF wneo5l
 


      RETURN
      END  SUBROUTINE  CH_neo




      SUBROUTINE CH_neo_single_grid_point(jm,pertrb)                             
! -------------------------------------------------------------------------- 
! -- get CHang_Hinton neoclassical diffusivites and related
! -- calls Ch_neo which adds to  dcoefp
! --------------------------------------------------------------------------

        USE nrtype,                                          ONLY : DP,I4B

        USE common_constants,                                ONLY : izero,zeroc,joupkev,                  &
                                                                    Electron_Rest_Mass,Electron_Charge,   &
                                                                    sqrt2,rootpi,root2pi,Proton_Mass,     &
                                                                    rootpio2

        USE dep_var,                                         ONLY : te,ti,angrot,en,etor,ene,             &
                                                                    rbp,dp4

        USE curden_terms,                                    ONLY : bp,q,dqdr,curden,shearp,eta,jboot

        USE grid_class,                                      ONLY : nj,r,dr,rcap,r2capi,xhm20,xhm2,       &
                                                                    fcap,gcap,hcap,r2cap,                 &
                                                                    xi11,xi110,xi33,xi330,ra,             &
                                                                    xips,xips0,r2cap,ravg_r,eps

        USE solcon_gcnmp,                                    ONLY : ntot,dudr,time,eqtime,                &
                                                                    include_neocl

        USE ions_gcnmp,                                      ONLY : nion,nprim,z,zeff,atw,dmassden,       &
                                                                    fi_index,nimp,zsq,dzdte

        USE neutral_beams,                                   ONLY : nbion

        USE plasma_properties,                               ONLY : mhd_dat,dischg

        USE MPI_data,                                        ONLY : myid !temporary

        IMPLICIT NONE

        INTEGER(I4B),INTENT(IN) :: jm,pertrb
        INTEGER(I4B) k,m,i,j
        REAL(DP) te_min,te_max !temp
        j = jm                              
        one_hundredth = 0.01_DP
        iwangrot = -3 ; ifus = 1
        xmasse  = Electron_Rest_Mass            ! kg
        xmassp  = Proton_Mass                   ! kg
        charge  = -Electron_Charge              ! ABS(qe) coul
        cesu    = charge*2.9979246e+09          ! esu
        charg4  = cesu**4                       ! esu**4
        enea    = 0.5 * ( ene(j)+ ene(j+1))     ! #/m**3
        tea     = 0.5 * (  te(j)+  te(j+1))     ! KEV
        tia     = 0.5 * (  ti(j)+  ti(j+1))
        curdna  = 0.5 * (curden(j)+curden(j+1)) ! amps/m**2
        dte     = (te(j+1)- te(j))/dr(j)
        dti     = (ti(j+1)- ti(j))/dr(j)
        zeffa   = 0.5_DP * (zeff(j)+zeff(j+1))
        zin     = 1.0_DP / zeffa
        angrota = 0.5_DP * (angrot(j)+angrot(j+1))  ! rad/sec
        enasum = zeroc
        DO k=1,nion
           za    (k)  = 0.5_DP * ( z(j,k)+ z(j+1,k))
           zsqa  (k)  = 0.5_DP * (zsq(j,k)+zsq(j+1,k))
           rzsqa (k)  = SQRT (zsqa(k))
           dzdtea(k)  = 0.5_DP * (dzdte(j,k)+dzdte(j+1,k))
           ena(k)     = 0.5_DP * (en(j,k)+en(j+1,k))
           enasum = enasum+ena(k)                                ! includes all ions
!              den(k)     = dudr(k,j)                             ! dudr was set in set_new_vars
           etaim(j,k) = zeroc
           dentia     = tia*(ene(j+1)-ene(j))/dr(j)
           IF (dentia .NE. zeroc  .AND. k .LE. nprim) &
                etaim(j,k) = enea * dti / dentia
        END DO


!        vionzgrd(j) = (angrot(j+1)*r2capi(j+1)/rcap(j+1) &
!             -angrot(j)*r2capi(j)/rcap(j))/dr(j)                 ! done in 
        !
        !  obtain delta and flux surface integrals
        !

        epsa   = (eps(j+1)+eps(j))*0.5_DP
        delta  = epsa
        rdelta = SQRT (delta)
        delta2 = delta**2
 
        CALL neointrp (xhm20, xhm2, j, xhm2a,time,eqtime)  
        CALL neointrp (xi110, xi11, j, xia(1,1),time,eqtime)
        CALL neointrp (xi330, xi33, j, xia(3,3),time,eqtime)
        CALL neointrp (xips0, xips, j, xipsa,time,eqtime)   

 

        xia(1,2) = xia(1,1)
        xia(1,3) = xia(3,3) / xhm2a
        xia(2,2) = xia(1,1)
        xia(2,3) = xia(1,3)

        !
        ! --- calculate trapped/circulating fraction HSJ
        !


        !
        ! --- see for example Hirshman et al N.F. 17,3 (1977) pg 611, eq. 8, for
        ! --- the following approximate form of xft
        xft = 1.0 - (1.0 - delta)**2 / SQRT (1.0 - delta2)/ &
             (1.0 + 1.46*rdelta)

        xfc     = 1.0 - xft
        ftfc(j) = 0.0
        IF (xfc .GT. 0.0)  ftfc(j) = xft/xfc
        !
        !  obtain geometric factors
        !
        del32  = delta*rdelta
        del32i = 1.0 / del32
        del52  = delta*del32
        gcapa  = 0.5 * ( gcap(j)+ gcap(j+1))
        fcapa  = 0.5 * ( fcap(j)+ fcap(j+1))
        hcapa  = 0.5 * ( hcap(j)+ hcap(j+1))
        r2capa = 0.5 * ( r2cap(j)+r2cap(j+1))
        grsq   = gcapa/r2capa

        !
        !  calculate bp, safac and shear PARAMETER
        !
        bp1       = rbp(j)/((r(j)+puny)* fcap(j)* gcap(j)* hcap(j))
        bp2       = rbp(j+1)/(r(j+1)* fcap(j+1)* gcap(j+1)* hcap(j+1))
        bpa       = 0.5 * (bp1+bp2)
        rbpa      = ra(j)*bpa
        safac     = ra(j)*mhd_dat%btor/(dischg%rmajor*bpa)
        safac     = ABS (safac)
!**        shearp(j) = ra(j) * ABS (q(j+1)-q(j))/(safac*dr(j))
        shearp(j) = ra(j)* dqdr(j)/safac
        shearp(j) = MAX (shearp(j), one_hundredth)

        !
        !     calculate Larmor radii and various collision frequencies
        !
        tejou   =  joupkev * tea                          ! joules
        teaerg  =  joupkev*1.e7_DP  * tea                 ! erg
        vthe    =  SQRT (2.0*tejou/xmasse)                ! m/sec
        rhoe    =  xmasse*vthe/(charge*bpa)               ! m 

        !     Coul Logarithm formula requires enea in 1/cm**3,te in ev:
        xlam    =  24.0_DP - DLOG (SQRT (enea*1.e-6)/(1.0D3*tea))
        xnue    =  1.33333e-6*ROOT2PI*enea*zeffa*charg4*xlam &
             / (SQRT (xmasse*1.e3*(teaerg)**3))    !1./sec
        xnuse(j) = sqrt2*dischg%rmajor*safac*xnue/(del32*vthe)   !coll freq /bounce freq
        tijou    = joupkev*tia
        tiaerg   = joupkev*1.e7_DP*tia
        xnubyv   = 1.33333_DP*ROOTPIO2*enea*1.e-6_DP*zeffa**3*charg4*xlam/tiaerg**2
        xnubyv   = xnubyv * 100._DP ! = 1./(vth tau) , 1./m
        xnusi(j) = sqrt2*dischg%rmajor*safac*xnubyv/del32        !coll freq /bounce freq

        DO  i=1,nion
           xmassi = atw(i)*xmassp
           vth(i) = SQRT (2.0_DP * tijou/ xmassi)                    !m/sec
           rho_ion(i) = xmassi*vth(i)/(rzsqa(i)*charge*bpa)           !m
           xnu(i) = 1.33333e-6*rootpi*ena(i)*zsqa(i)**2*charg4*xlam &
                / SQRT (xmassi*1.e3*tiaerg**3)                   !1/sec
           xnus(i,j) = sqrt2*dischg%rmajor*safac*xnu(i)/(del32*vth(i)) 
        ENDDO
       !
       !  calculate fit parameters of Rawls, Chu, and Hinton for given zeff
       !  NEED to check  units on this:
        DO m=1,3
           DO k=m,3
              xk0(m,k) = f1i(zin,xk01(m,k),xk0i(m,k))
              a(m,k)   = f14(zin,a1(m,k),a4(m,k))
              b(m,k)   = f14(zin,b1(m,k),b4(m,k))
              c(m,k)   = f14(zin,c1(m,k),c4(m,k))
              IF (k .NE. 3)THEN
                 xk(m,k) = xk0(m,k) &
                      * ((1.0/(1.0 + a(m,k) * SQRT (xnuse(j))+b(m,k)*xnuse(j))) &
                      * xia(m,k)/(1.95 * rdelta) + 0.5 * c(m,k)**2*xnuse(j)*delta*xipsa &
                      / (b(m,k)*(1.0 + c(m,k)*xnuse(j)*del32)))
              ELSE
                 xk(m,k) = &
                      xk0(m,k)*(1.0/(1.0 + a(m,k) * SQRT (xnuse(j))+b(m,k)*xnuse(j))) &
                      *(1.0 / (1.0 + c(m,k)*xnuse(j)*del32))*xia(m,k) / (1.95 * rdelta)
              ENDIF
           ENDDO
        ENDDO
        !
        !  calculate resistivity
        !  note that for zeffa >3.0 the conductivity reduction
        !  due to e-e collisons is assumed negligible.
        !  neoresist returns effective trapped particle fraction
        !  ftrap,etaspitzer,eta(j),and xkap11. (xkap11 is numerically
        !  essentially identical to the xkap11 defined above.)
        !
        te_min =Minval(te)
        te_max =Maxval(te)
        IF(te_min .lt. 0.00)THEN
           write(myid+100,FMT='("neo_tport reports: te min lt 0.0")')
           write(myid+100,8)te
8        format(5(2x,1pe14.6))
           flush(myid+100)
        ENDIF
        IF(te_max .gt. 1.e2)THEN
           write(myid+100,FMT='("neo_tport reports: te max gt 100")')
           write(myid+100,9)te
9        format(5(2x,1pe14.6))
           flush(myid+100)
        ENDIF

        CALL neoresist (zeffa,wneot,xft, &
             xnue,xnuse(j),xmasse,enea,charge, &
             ftrap(j),etaspitzer,eta(j),xkap11)
        eta(j) = wneo(4,4)*eta(j)   !ohm m
        sumdz = 0.0
        DO  k=1,nion
           sumdz = sumdz+dzdtea(k)*ena(k)
        ENDDO
        sumdz = sumdz*tea/enea

        !
        ! ----------------------------------------------------------------------
        ! ---                    neoclassical transport                      ---
        ! --- Note that if background momentum diffusivity is taken as neoclassical
        ! --- ion energy diffusivity, chiineo. Hence we may have to call
        ! --- CH_neo to get chiineo even if neoclassical is not used otherwise.
        ! --- Ch_neo loads  dcoefp 
        ! ----------------------------------------------------------------------
        !
  
        IF (include_neocl  .EQ. 1 .OR. iwangrot .EQ. -3)                          &
             CALL CH_neo (j,sumdz,xnuse,xnusi,xnus,ftrap,xk,xhm2a,ena,enasum,     & 
             za,zeffa,wneo,fcapa,hcapa,rbpa,wneo1,wneo2,wneo3,wneo4,wneo5,        & 
             dwifs,tia,tea,teaerg,enea,del32,delta,rdelta,xipsa,bpa,xnue,         & 
             rhoe,rho_ion,xkap11,xnu,ifus,zsqa,xia,jboot,iwangrot,delta2,del32i)    
                                                        
        RETURN

      END SUBROUTINE Ch_neo_single_grid_point




      SUBROUTINE neoresist (zeffa, wneot, xft, xnue, xnuse, xmasse,     &
                            enea, charge, fftrap, etaspitzer,           &
                            etaparallel, xkap11)
!
! ----------------------------------------------------------------------
! subroutine calculates the neoclassical resistivity with finite
! inverse aspect ratio and collisional effects.
! see Hirshman et al.,N.F. 17,3,(1977)pg. 611
!                     N.F. 21,9,(1981)pg. 1158,eq. 7.41
!  Note that this formulation is valid only for zeff< 3.0
!  for zeff > = 3.0 we set cr (the conductivity reduction due to
!  e-e collisions) equal to 0.0
!
! --- input
!
!  zeffa        zeff on half grid
!  xft          trapped particle fraction (collisionless)
!  wneot        multiplier,normally set to 1.0,for effective trapped
!               electron fraction.
!  xnue         e-e collision frequency (includes zeff term)
!  xnuse        collision/bounce frequency (includes zeff)
!  xmasse       electron mass,kgrams
!  enea         electron density, m-3, on half grid
!  charge       electron charge, coul
!
! --- output
!  fftrap        effective trapped electron fraction
!
!  etaparallel   parallel resistivity (in ohm m)
!
!  etaspitzer    Spitzer resistivity (in ohm m)
!  xkap11        zeff dependent factor of etaspitzer (note that zeff
!                also appears in defn of xnue,xnuse)
! ------------------------------------------------------------------ HSJ
!
      USE nrtype
      USE error_handler,           ONLY : iomaxerr, lerrno,terminate
      USE io_gcnmp,                ONLY : nlog
      IMPLICIT NONE
!
      REAL(DP) xft,zeffa,fftrap,etaspitzer,etaparallel,xnue,xnuse,  &
             xnusem,enea,charge,xmasse,z,cr,zeta,zr,ftad,zta,cra,   &
             wneot,xkap11,xmasse_g,enea_cm3,charge_esu

!
      cr(z)   = 0.56*(3.0-z)/(3.0+z)/z
      zeta(z) = 0.58+0.20*z


!
!  --- eq. 7.36, hirshman N.F. 1981, pg 1157
!
      zr(z) = ((0.222*z+1.198)*z+1.0)/((0.753*z+2.966)*z+1.0)
!
!     note: tauee = zeffa/xnue
!


      xmasse_g = xmasse*1000.
      enea_cm3 = enea*1.e-6
      charge_esu = charge*2.9979246e+09

      xnusem     = xnuse/zeffa
      xkap11     = 1.0 / zr(zeffa)
      etaspitzer = xmasse_g*xnue/(enea_cm3*charge_esu**2)         ! in sec
      etaspitzer = etaspitzer/xkap11 ! eq 7.36, pg 1157, NF 1981
      zta        = zeta(zeffa)
      IF (zeffa .LT. 3.0) THEN
        cra = cr(zeffa)
      ELSE
        cra = 0.0
      END IF
      ftad        = xft/(1.0 + zta*xnusem)
      fftrap      = wneot*ftad*(1.0 + cra*(1.0-ftad))
      etaparallel = etaspitzer/(1.0-fftrap)
      IF(etaparallel .LT. 0.0)THEN
         Print *,'xft,xnusem,wneot=',xft,xnusem,wneot
         Print *,'xnue=',xnue
         PRINT *,'sub neoresist, fftrap , etaspitzer =',  fftrap,etaspitzer
         lerrno = iomaxerr + 53
         CALL terminate (lerrno,nlog)
      ENDIF
      !convert to ohm m 
      etaparallel = etaparallel * 8.98755179e+9   ! ohm m
      etaspitzer  = etaspitzer  * 8.98755179e+9   ! ohm m
      RETURN
      END SUBROUTINE neoresist



  SUBROUTINE neointrp (x0, x1, j, xav,time,eqtime)
!
      USE nrtype,          ONLY : I4B,DP

      USE common_constants,       ONLY : zeroc

      IMPLICIT  NONE 

      INTEGER(I4B) j 

      REAL(DP)   x0(*), x1(*),xav,x0av,x1av,time,dteq,eqtime,eqtim0

      dteq = zeroc    ; eqtim0 = zeroc     !no time dep mhd  
      x0av = 0.5 * (x0(j) + x0(j+1))
      xav  = x0av
      IF (dteq .GT. zeroc ) THEN
        x1av = 0.5 * (x1(j) + x1(j+1))
        xav  = x0av + (time-eqtim0) * (x1av-x0av) / dteq
      END IF
      RETURN
  END SUBROUTINE neointrp

 
  SUBROUTINE neo_allocate
     USE grid_class,                             ONLY : nj
     USE ions_gcnmp,                             ONLY : nion
     USE common_constants,                       ONLY : zeroc
     IMPLICIT NONE
     
          IF(.NOT. ALLOCATED(rho_ion))ALLOCATE(rho_ion(nion)) 
          IF(.NOT. ALLOCATED(xnu))ALLOCATE(xnu(nion)) 
          IF(.NOT. ALLOCATED(xnus))ALLOCATE(xnus(nion,nj)) 
          IF(.NOT. ALLOCATED(xnuse))ALLOCATE(xnuse(nj)) 
          IF(.NOT. ALLOCATED(ftrap))ALLOCATE(ftrap(nj))  
!          IF(.NOT. ALLOCATED(slpres))ALLOCATE(slpres(nj)) 
          IF(.NOT. ALLOCATED(za))ALLOCATE(za(nion))
          IF(.NOT. ALLOCATED(zsqa))ALLOCATE(zsqa(nion))
          IF(.NOT. ALLOCATED(rzsqa))ALLOCATE(rzsqa(nion))
          IF(.NOT. ALLOCATED(dzdtea))ALLOCATE(dzdtea(nion))
          IF(.NOT. ALLOCATED(ena))ALLOCATE(ena(nion))
!          IF(.NOT. ALLOCATED(den))ALLOCATE(den(nion))
          IF(.NOT. ALLOCATED(etaim))ALLOCATE(etaim(nj-1,nion))
          IF(.NOT. ALLOCATED(vionzgrd))ALLOCATE(vionzgrd(nj))
          IF(.NOT. ALLOCATED(ftncl))ALLOCATE(ftncl(nj))
          IF(.NOT. ALLOCATED(ftfc))ALLOCATE(ftfc(nj))
          IF(.NOT. ALLOCATED(xnusi))ALLOCATE(xnusi(nj))
          IF(.NOT. ALLOCATED(vth))ALLOCATE(vth(nion))
          IF(.NOT. ALLOCATED(dwifs))ALLOCATE(dwifs(nj))

          xnus(:,:)    = zeroc         ; xnu(:)      = zeroc
          xnuse(:)     = zeroc         ; ftrap(:)    = zeroc
          etaim(:,:)   = zeroc         ; vionzgrd(:) = zeroc
          ftfc(:)      = zeroc         ; xnusi(:)    = zeroc
          ftncl(:)     = zeroc         ; dwifs(:)    = zeroc
          vth(:)       = zeroc         ; vth(:)      = zeroc
     RETURN
   
  END   SUBROUTINE neo_allocate



   END MODULE neo_tport  
