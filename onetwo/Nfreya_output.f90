

  SUBROUTINE Nfreya_output(timet, t)
!----------------------------------------------------------------------
! -- This routine is used for Nfreya and P_nfreya
! -- USES external routines  header,trapv
! -- writes to units nqik,nout
! -- INPUT
!       timet
!       t
!       qbbe(:,:,:) 
!       qbbi(:,:,:) 
!
!---------------------------------------------------------------HSJ----
    USE nrtype,                                   ONLY : DP,I4B

    USE param,                                    ONLY : ke,kb,kj

    USE numbrs,                                   ONLY : nj

    USE Plasma_properties,                        ONLY : neut_beam
  
    USE common_constants,                         ONLY : zeroc,izero

    USE geom,                                     ONLY : hcap,volfac

    USE mesh,                                     ONLY : r,roa

    USE io,                                       ONLY : jprt,nout,nqik

    USE nub,                                      ONLY : pbeam,fap,fbcur,  &
                                                         fwall,ebeam,qbbe, &
                                                         qbbi,qb,ibeam,    &
                                                         hibr,zeta,nbeams
    USE P_Nfreya_12_interface,                    ONLY : use_P_Nfreya

    USE nub2,                                     ONLY : bion,bneut,forb,    &
                                                         qbsav,taupb,taueb,  &
                                                         fber,fb00,fb01,fb10,&
                                                         fb11,wb00,wb01,wb11,&
                                                         wb10,hdep,fbe,fbi

    USE transp,                                   ONLY : beam_data 

    USE yoka 

    IMPLICIT NONE
    

    INTEGER(I4b) j,jb,j1prt,ic

    REAL(DP) pbap,pbsap,pwall,pborb,pbplaf,pbplas,pbel,pbion
    REAL(DP) fsap,fw,florb,fpbcx,fpbe,fpbi,fxorb
    REAL(DP) timet,t
    REAL(DP) pbeamf(ke,kb),pbeams(ke,kb) 
    REAL(DP) fpe(ke,kb),fpi(ke,kb),fpcx(ke,kb)
    REAL(DP) psum(kj), pdum(kj), qdum(kj)
    LOGICAL  nf_printout 

      pbap   = zeroc
      pbsap  = zeroc
      pwall  = zeroc
      pborb  = zeroc
      pbplaf = zeroc
      pbplas = zeroc
      pbel   = zeroc
      pbion  = zeroc
      pfil   = zeroc
      sthru  = zeroc
      nf_printout = .FALSE.

      IF(use_P_Nfreya)THEN ! check units
         nbeams = neut_beam%nbeams 
 !do jb=1,nbeams
 !  do ic =1,ke 
 !      print *,'nb%fap=',neut_beam%fap(ic,jb),ic,jb
 !  enddo
 !enddo
         Do jb =1,nbeams 

            IF(beam_data%beamlet_active(jb) .NE. izero)THEN
               nf_printout = .TRUE.
               pbeam(1:ke,jb)   = neut_beam%pbeam(1:ke,jb)  
               fap(1:ke,jb)     = neut_beam%fap(1:ke,jb)
               fwall(1:ke,jb)   = neut_beam%fwall(1:ke,jb)
               forb(1:ke,jb)    = neut_beam%forb(1:ke,jb)
               ebeam(1:ke,jb)   = neut_beam%ebeam(1:ke,jb)
               bneut(1:ke,jb)   = neut_beam%bneut(1:ke,jb)
               bion(1:ke,jb)    = neut_beam%bion(1:ke,jb)
               hibr(1:nj,1:ke,jb)  = neut_beam%hibr(1:nj,1:ke,jb)
               hdep(1:nj,1:ke,jb)  = neut_beam%hdep(1:nj,1:ke,jb) 
               zeta(1:nj,1:ke,jb)  = neut_beam%zeta(1:nj,1:ke,jb)
               fb11(1:ke,jb)    = neut_beam%fb11(1:ke,jb)
               fb10(1:ke,jb)    = neut_beam%fb10(1:ke,jb)
               fb01(1:ke,jb)    = neut_beam%fb01(1:ke,jb)
               fb00(1:ke,jb)    = neut_beam%fb00(1:ke,jb)
               wb11(1:ke,jb)    = 100._DP*neut_beam%wb11(1:ke,jb)    ! units of cm 
               wb10(1:ke,jb)    = 100._DP*neut_beam%wb10(1:ke,jb)
               wb01(1:ke,jb)    = 100._DP*neut_beam%wb01(1:ke,jb)
               wb00(1:ke,jb)    = 100._DP*neut_beam%wb00(1:ke,jb)
               fber(1:ke,jb)    = neut_beam%fber(1:ke,jb)
            ELSE 
               pbeam(1:ke,jb)   = zeroc
               fap(1:ke,jb)     = zeroc
               forb(1:ke,jb)    = zeroc 
               ebeam(1:ke,jb)   = zeroc 
               bneut(1:ke,jb)   = zeroc
               bion(1:ke,jb)    = zeroc
               hibr(:,1:ke,jb)  = zeroc
               hdep(:,1:ke,jb)  = zeroc
               zeta(:,1:ke,jb)  = zeroc
               fb11(1:ke,jb)    = zeroc
               fb10(1:ke,jb)    = zeroc
               fb00(1:ke,jb)    = zeroc
               fb01(1:ke,jb)    = zeroc
               wb11(1:ke,jb)    = zeroc
               wb10(1:ke,jb)    = zeroc
               wb00(1:ke,jb)    = zeroc
               wb01(1:ke,jb)    = zeroc
               fber(1:ke,jb)    = zeroc
            ENDIF
         ENDDO
      ENDIF




      IF(.NOT. use_P_Nfreya .AND. ibeam .LT. 2 ) RETURN 

      IF( use_P_Nfreya .AND. .NOT. nf_printout ) RETURN 



      DO  jb=1,nbeams
          IF(use_P_Nfreya)THEN
             ! beamlet_active not defined if use_P_nfreya = False
             IF(beam_data%beamlet_active(jb) .EQ. izero)CYCLE 
          ENDIF
!         IF(use_P_Nfreya .AND. beam_data%beamlet_active(jb) .EQ. izero)CYCLE ! cant use in ifort
         DO  ic=1,3
            pbap   = pbap  + pbeam(ic,jb)
            pbsap  = pbsap + fap(ic,jb)*pbeam(ic,jb)
            pwall  = pwall + fwall(ic,jb)*pbeam(ic,jb)
            pborb  = pborb + forb(ic,jb)*pbeam(ic,jb)
            pbeamf(ic,jb) = pbeam(ic,jb)* &
                        (1.0_DP-fap(ic,jb)-fwall(ic,jb)-forb(ic,jb))
            pbplaf = pbplaf + pbeamf(ic,jb)
            fpe(ic,jb)  = zeroc
            fpi(ic,jb)  = zeroc
            fpcx(ic,jb) = zeroc
            DO  j=1,nj
               psum(j) = qbbe(j,ic,jb)
               pdum(j) = qbbi(j,ic,jb)
            ENDDO
            CALL trapv (r,psum,hcap,nj,fpe(ic,jb))
            CALL trapv (r,pdum,hcap,nj,fpi(ic,jb))
            CALL trapv (r,qb(1,ic,jb),hcap,nj,pbeams(ic,jb))
            pbeams(ic,jb) = volfac*pbeams(ic,jb)
            IF (pbeams(ic,jb) .GT. zeroc)THEN
               fpe(ic,jb)  = volfac*fpe(ic,jb)/pbeams(ic,jb)
               fpi(ic,jb)  = volfac*fpi(ic,jb)/pbeams(ic,jb)
               fpcx(ic,jb) = 1.0_DP-fpe(ic,jb)-fpi(ic,jb)
               pbplas      = pbplas + pbeams(ic,jb)
               pbel        = pbel   + fpe(ic,jb)*pbeams(ic,jb)
               pbion       = pbion  + fpi(ic,jb)*pbeams(ic,jb)
               pfil        = pfil   + fpcx(ic,jb)*pbeams(ic,jb)
            ENDIF
         ENDDO
      ENDDO
      fsap   = zeroc
      fw     = zeroc
      florb  = zeroc
      IF (pbap .NE. zeroc)THEN  
         fsap   = pbsap/pbap
         fw     = pwall/pbap
         florb  = pborb/pbap
      ENDIF

      fpbe   = zeroc

      fpbi   = zeroc
      fpbcx  = zeroc
      IF (pbplas .NE. zeroc)THEN
         fpbe   = pbel/pbplas
         fpbi   = pbion/pbplas
         fpbcx  = pfil/pbplas
      ENDIF

      ptor   = pbap - pbsap

      IF (ptor .GT. zeroc)  sthru = (pwall/ptor) * 100.0_DP
      pfil   = pfil * 1.0e-6

! print out neutral beam injection sources

      DO  jb=1,nbeams
          IF(use_P_Nfreya)THEN
             ! beamlet_active not defined if use_P_nfreya = False
             IF(beam_data%beamlet_active(jb) .EQ. izero)CYCLE 
          ENDIF
!         IF(use_P_Nfreya .AND. beam_data%beamlet_active(jb) .EQ. izero)CYCLE
         DO  ic=1,3

            CALL header (nout, timet, t)
            WRITE (nout, 1146)  jb, ic
            DO  j=1,nj
               j1prt = ((j-1)/jprt)*jprt
               IF (.NOT.((j1prt .NE. j-1) .AND. (j .NE. nj)))THEN

                  WRITE (nout,1153) j,r(j),hibr(j,ic,jb),hdep(j,ic,jb), &
                       zeta(j,ic,jb), qbsav(j,ic,jb), qb(j,ic,jb), &
                       fbe(j,ic,jb), fbi(j,ic,jb), taupb(j,ic,jb), taueb(j,ic,jb)
               ENDIF
            ENDDO

            WRITE (nout,8110)         ebeam(ic,jb)
            WRITE (nout,8120)         bion(ic,jb)
            WRITE (nout,8130)         bneut(ic,jb)
            WRITE (nout,8140) pbap,   pbeam(ic,jb)
            WRITE (nout,8150) fsap,   fap(ic,jb)
            WRITE (nout,8160) fw,     fwall(ic,jb)
            WRITE (nout,8161) florb,  forb(ic,jb)
            WRITE (nout,8162) pbplaf, pbeamf(ic,jb)
            WRITE (nout,8163) pbplas, pbeams(ic,jb)
            WRITE (nout,8164) fpbe,   fpe(ic,jb)
            WRITE (nout,8166) fpbi,   fpi(ic,jb)
            WRITE (nout,8168) fpbcx,  fpcx(ic,jb)
            WRITE (nout,8170)
 
            fxorb = forb(ic,jb)/(1.0_DP-fap(ic,jb)-fwall(ic,jb))

            WRITE (nout,8180) fb11(ic,jb), wb11(ic,jb), &
                   fb10(ic,jb), wb10(ic,jb), fb01(ic,jb), wb01(ic,jb), &
                   fb00(ic,jb), wb00(ic,jb), fxorb, fber(ic,jb)
         ENDDO
 
      ENDDO

      CALL header (nqik, timet, t)

      ! prints out all beams. For P_Nfreya beams calcs the source
      ! will be zero for beams that are off.
      WRITE (nqik,8200)
      WRITE (nqik,8110)        (( ebeam(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8120)        ((  bion(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8130)        (( bneut(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8140) pbap,  (( pbeam(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8150) fsap,  ((   fap(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8160) fw,    (( fwall(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8161) florb, ((  forb(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8162) pbplaf,((pbeamf(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8163) pbplas,((pbeams(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8164) fpbe,  ((   fpe(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8166) fpbi,  ((   fpi(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8168) fpbcx, ((  fpcx(ic,jb),ic=1,3),jb=1,nbeams)
      WRITE (nqik,8170)


      DO jb=1,nbeams
          IF(use_P_Nfreya)THEN
             ! beamlet_active not defined if use_P_nfreya = False
             IF(beam_data%beamlet_active(jb) .EQ. izero)CYCLE 
          ENDIF
!         IF(use_P_Nfreya .AND. beam_data%beamlet_active(jb) .EQ. izero)CYCLE
         WRITE (nqik,8210) jb,(ebeam(ic,jb),ic=1,3)
         DO  j=1,nj
            j1prt = ((j-1)/jprt)*jprt
            IF (.NOT.(j1prt .NE. j-1 .AND. j .NE. nj))THEN
               WRITE (nqik,8220) j,r(j),roa(j),    &
                    (qb(j,ic,jb),fbe(j,ic,jb),fbi(j,ic,jb),ic = 1,3)
            ENDIF
         ENDDO

      ENDDO

 8110 FORMAT (// ' particle energy (keV)', 22x, 'total',6(4x,f8.3))
 8120 FORMAT (' ion beam intensity (part./s)',24x,1p6e12.3)
 8130 FORMAT (' neutral beam intensity (part./s) to ap.',13x,1p6e12.3)
 8140 FORMAT (' neutral beam power (W) to aperture',6x,1p7e12.3)
 1146 FORMAT (20x,'neutral beam injection sources, beam number ',i2, &
            /30x,'component',i6, &
            //15x,2( '  normalized'),4x,'average', &
            5x,'fast ion',6x,'delayed',8x,'energy fraction', &
            5x,'p. slowing',3x,'e. slowing' &
            /4x,'j',7x,'r',1x,2(5x,'hot ion'), '  pitch angle', &
            2(3x,'e. source*'),8x,'deposited in', &
            4x,2(4x,'down time') &
            /10x,'(cm)',3x,'birth rate',3x,'dep. rate',5x,'cosine',1x, &
            2(4x,'(W/cm**3)'),4x,'electrons',6x,'ions', &
            2(8x,'(s)',2x))
 8150 FORMAT (' fraction stopped by aperture',9x,7(3x,f9.4))
 1153 FORMAT (1x,i4,f9.2,3f12.4,1x,1p2e13.3,0p2f12.4,1p2e13.3)
 8160 FORMAT (' fraction incident on wall (shinethrough)', f9.4,6(3x,f9.4))
 8161 FORMAT (' fraction lost on orbits',14x,7(3x,f9.4))
 8162 FORMAT (' neutral beam power (W) in plasma ',7x,1p7e12.3)
 8163 FORMAT (' slowed  beam power (W) in plasma*',7x,1p7e12.3)
 8164 FORMAT (' fraction deposited in electrons',6x,7(3x,f9.4))
 8166 FORMAT (' fraction deposited in ions',11x,7(3x,f9.4))
 8168 FORMAT (' fraction lost by fast ion charge ex.',1x, 7(3x,f9.4))
 8170 FORMAT (/' * excludes aperture, shinethrough, and orbit', &
             ' losses, but includes loss due to fast ion charge ex.')
 8180 FORMAT (/ ' fraction and width of various type orbits' // &
           15x, 'type',16x,'fraction',7x,'width'/51x,'(cm)' /   &
            3x, 'passing and axis-circling',3x,f12.4,f12.2  /   &
            3x, 'passing and not circling',4x,f12.4,f12.2   /   &
            3x, 'trapped and axis-circling',3x,f12.4,f12.2  /   &
            3x, 'trapped and not circling',4x,f12.4,f12.2   /   &
            3x, 'lost on orbit',15x,f12.4                   /   &
            3x, 'error detected',14x,f12.4)
 8200 FORMAT (   20x, 'neutral beam injection sources')
 8210 FORMAT (// 20x, 'fast ion power source and fraction deposited', &
                    ' in electrons and ions    (beam',i3,')' // &
                      3x,'energy (keV):',12x,3(f8.3,24x) / &
                      4x,'j',4x,'r',3x,'r/a',5x,&
                      3(1x,'W/cm**3',5x,'elec',4x,'ions',7x))
 8220 FORMAT (1x,i4,f6.1,f5.1,3(2x,1pe12.3,0p2f8.3,2x))
      RETURN

  END SUBROUTINE Nfreya_output
