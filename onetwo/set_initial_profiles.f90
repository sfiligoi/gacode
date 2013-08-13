      SUBROUTINE set_initial_profiles
! ----------------------------------------------------------------------
!    generate initial  profiles for all dependent variables.
!    The statefile option can supply initial profiles or
!    they can be obtained from inone
!    IF the statefile option is selected time dependent profiles
!    beyond the inital time  are input in inone.
!    Note that the inone file values for the initial time are
!    overwritten by the statefile values if the statefile option is selected.
! ----------------------------------------------------------------------
!
      USE nrtype,                         ONLY : I4b,DP

      USE param,                          ONLY : ksplin,kprim,kimp,kbctim,krf,kj

      USE iterdbmd,                       ONLY : initialize_from_statefile,         &
                                                 statefile_name

      USE bd_condtn,                      ONLY : rnormin,njte,njti,njene,njzef,     &
                                                 renpin,reniin,renein,enein,        &
                                                 iprofnbr,eni,enp,enpb,bctime,      &
                                                 knotsene,gamene,alpene,            &
                                                 enec,eneb,gameni,alpeni,           &
                                                 enic,enib,enpc,enpb,alpenp,gamenp, &
                                                 tec,teb,alpte,gamte,               &
                                                 tic,tib,alpti,gamti,               &
                                                 xjc,xjb,alpxj,gamxj,               &
                                                 zeffc,zeffb,alpzef,gamzef,         &
                                                 rtein,enein,tein,knotste,          &
                                                 rtiin,tiin,knotsti,mtimeprf,       &
                                                 rcurdein,curdenin,knotscur,        &
                                                 rzeffin,zeffin,knotszef,           &
                                                 knotscurb_external, rcurb_external,&
                                                 njcurb_external,curbeam_external,  &
                                                 bparenp, bpareni,bparte, bparti,   &
                                                 bparang, bparene,bparzeff,         &
                                                 bparcur,bparkpol,                  &
                                                 bparcurb_external,taupin,          &
                                                 en_bc_inpt,ren_bc_inpt,            &
                                                 ken_bc_inpt,density_mult,          &
                                                 zeff_mult,ene_mult,te_mult,        &
                                                 ti_mult, ang_mult

      USE tordlrot,                       ONLY : rangrot,angrotin,angrot,           &
                                                 njinang,iangrot,knotsang
                                                 
      USE numbrs,                         ONLY : nj,nprim,nimp,nion,nbctim

      USE ions,                           ONLY : njenp,njeni,zeff,z,zsq,namei

      USE io,                             ONLY : ncrt,nout

      USE mesh,                           ONLY : r,roa

      USE soln,                           ONLY : te,ti,ene,en,curden,etor,&
                                                 rbp,diffeq_methd,njcur

      USE solcon,                         ONLY : time,external_beam_cur

      USE sourc,                          ONLY : njqin,qine,qini,curb_external

      USE common_constants,               ONLY : zeroc

      USE rf,                             ONLY : rfmode

      USE string_util,                    ONLY : to_upper_case1

      USE flags,                          ONLY : inenez,zfrac

      USE nub2,                           ONLY : enbeam

      USE fusion,                         ONLY : enalp,dtnfus

      USE gpsi,                           ONLY : pppsi_eqdsk,presspsi_eqdsk, &
                                                 fpsi_eqdsk,ffppsi_eqdsk,    &
                                                 qpsi_eqdsk,nxeqd_eqdsk

      USE psig,                           ONLY : fpsi,ffppsi,presspsi,qpsi,pppsi


      IMPLICIT NONE

      INCLUDE 'imsl.i'


      INTEGER(I4b) i,jnormin,j,k,njin,m,ierr,ierr1,npsig
      REAL(DP) pos_def
!      REAL(DP),dimension(:)       ::   profin(nj),xdum(nj),ydum(nj)
      REAL(DP) profin(kj),xdum(kj),ydum(kj)
      CHARACTER*8           tstring
      LOGICAL ene_set,te_set,ti_set,zeff_set,enp_set,eni_set, &
              curden_set,angrot_set
      pos_def = 5._DP*EPSILON(0.0_DP)
      imslmd = 'init- sp'
      m          = 1

      ene_set    = .FALSE.
      te_set     = .FALSE.
      ti_set     = .FALSE.
      zeff_set   = .FALSE.
      enp_set    = .FALSE.
      eni_set    = .FALSE.
      curden_set = .FALSE.
      angrot_set = .FALSE.
      ierr       = 0


      mtimeprf = 1
      do i=2,nbctim
        mtimeprf  = i
        if (bctime(i) .ge. time) go to 2
      end do

2      IF(initialize_from_statefile) THEN

         CALL read_statefile  ! read   text or netcdf file
                              ! and copy  to onetwo variables
                              ! (in sub set_onetwo_vars)

         IF(nj .GT. kj)THEN
            ! nj read from statefile is incompatible with this onetwo max array size of kj
            Write(nout,5)statefile_name(1:LEN_TRIM(statefile_name)),nj,kj
            Write(ncrt,5)statefile_name(1:LEN_TRIM(statefile_name)),nj,kj
5           FORMAT(2x," ERROR in using statefile ",a,/,  &
                   2x," size of arrays in statefile = ",i5,/, &
                   2x," Onetwo array size kj =",i5,/,        &
                   2x," kj must be >= nj, (recompile code with bigger kj)")
            ierr = 1
         ENDIF
         IF(ierr == 1) CALL EXIT

         !set primary ions  consistent with statefile condition:
         ! Note that since primary and impurity ions do not have
         ! time dependent input densities they are handled differently
         ! than the other inputs below.
               DO   i=1,nprim
                  renpin(1:nj,i) =  roa(1:nj)  !note i is species not time  here
                  njenp(i)       =  nj
                  enp(1:nj,i)    =  en(1:nj,i) 
               ENDDO
               enp_set = .TRUE.
         !set impurities  consistent with statefile condition:
               DO i=1,nimp
                  k = nprim + i
                  reniin(1:nj,i) =  roa(1:nj)  !note i is species not time  here
                  njeni(i)       =  nj
                  eni(1:nj,i)    =  en(1:nj,k) 
               ENDDO
               eni_set = .TRUE.
         ! remaining profiles are zeroed in the INITIAL time slot
         ! because they may have been read in in inone but are to be replaced
         ! by statefile input. For example, tein(1:nj,1) is the profile
         ! for te for the initial time. This profile will be supplied
         ! by the statefile below.



         !set ene  consistent with statefile condition
              renein(1:nj,1)   =  roa(1:nj)
              njene            =  nj
              enein(1:nj,1)    =  ene(1:nj)
              knotsene(1)      =  nj
              DO j=2,nbctim
                 IF(knotsene(j) .lt. 3)THEN
                    knotsene(j) = knotsene(1)
                    renein(1:nj,j)   = renein(1:nj,j-1)
                    enein(1:nj,j)    = enein(1:nj,j-1)
                 ENDIF
              ENDDO
              eneb(1)          =  ene(nj)
              ene_set          =  .True.
         !set TE consistent with statefile condition:
              rtein(1:nj,m)    = roa(1:nj)
              njte             = nj
              tein(1:nj,1)     = te(1:nj)
              knotste(1)       = nj
              te_set           = .TRUE.
         !set TI consistent with statefile condition:
              rtiin(1:nj,m)    = roa(1:nj)
              njti             =  nj
              tiin(1:nj,1)     =  ti(1:nj)
              knotsti(1)       =  nj
              ti_set           = .TRUE.
         !set curden  consistent with statefile condition:
              rcurdein(1:nj,m) = roa(1:nj)
              njcur            = nj
              curdenin(1:nj,1) = curden(1:nj) 
              knotscur(1)      = nj
              curden_set       = .TRUE.
         !set angrot  consistent with statefile condition:
              rangrot(1:nj,m)  = roa(1:nj)
              njinang          = nj
              angrotin(1:nj,1)    = angrot(1:nj)
              knotsang(1)      = nj
              angrot_set       = .TRUE.
         !set zeff  consistent with statefile condition:
              rzeffin(1:nj,m)  = roa(1:nj)
              njzef            = nj
              zeffin(1:nj,1)   = zeff(1:nj)  
              knotszef(1)      = nj
              zeff_set         = .TRUE.

         ! NOTE: eqdsk type quantities qpsi,ffppsi,pppsi,psival,presspsi,fpsi
         ! were set to
         ! pppsi_eqdsk,presspsi_eqdsk,fpsi_eqdsk,ffppsi_eqdsk,qpsi_eqdsk
         ! which were loaded from mhd_dat%** state file input in set_12_gcnmp_vars
         ! units and order of values were converted there
         ! here we move them into the local eqdsk values used in onetwo.
 
            !IF( .NOT. ALLOCATED(presspsi_eqdsk)) &
            !    CALL STOP(" sub set_initial profiles, error in eqdsk value init",1)
            !
            !     fpsi(:)     = fpsi_eqdsk(:)  
            !     ffppsi(:)   = ffppsi_eqdsk(:)
            !     presspsi(:) = presspsi_eqdsk(:)
            !     qpsi(:)     = qpsi_eqdsk
            !     pppsi(:)    = pppsi_eqdsk
            !     psival(1:npsi) = psival_eqdsk(:)

       ENDIF  ! statefile  initialization


!     Check for profiles WITH input DATA to spline
!     check that rnormin array can be used IF necessary (jnormin .GE. 3):
!
      DO j=1,ksplin
        jnormin = j
        IF (ABS(rnormin(j)- 1.0_DP) .LT.  pos_def )  go to 1
      END DO
      jnormin = 0
1     IF (jnormin .LE. 3)  jnormin = 0
      IF (jnormin .GT. kj) jnormin = 0



      !check for parabolic input by sensing on njin.
      !parabolic input is used if njin remains 0:
      njin = 0
      DO  i=1,nprim
         IF (njenp(i) .GT. 0) njin = njenp(i)
      ENDDO
      DO   i=1,nimp
         IF (njeni(i) .GT. 0)  njin = njeni(i)
      ENDDO






!----------------------------------------------------------------------
! --- initial electron density: 
! --- inenez = +1 or -1 means ene and zeff are input, enp, eni are calcualted
! --- inenez = 0 means determine ene zeff from input of primaries and impurities
! --- inenez = -98,-99 see below
!---------------------------------------------------------------------
    IF (njene .GT. 0 .AND. inenez .NE.  0 )THEN
      ierr1 = 1
      IF (ABS(renein(2,m)) .LT. pos_def .AND. (jnormin .GE. 3) .AND. ene_set == .FALSE. ) THEN
        renein(1:jnormin,m) = rnormin(1:jnormin)
        njene = jnormin
        ierr1 = 0
      ELSE IF (renein(2,m) .EQ. 0.0 .AND. (jnormin .LT. 3) .AND. ene_set == .FALSE.) THEN
        ierr1 = 1
      ELSE
        CALL knotnum(renein(1,m),njene,ksplin)
        IF (njene .GE. 3)  ierr1 = 0
      END IF
      IF (ierr1 .EQ. 1) THEN
        WRITE (ncrt, 40) 'ne'
        WRITE (nout, 40) 'ne'
      END IF
      ierr     = MAX0 (ierr, ierr1)
      iprofnbr = kprim + kimp + 1
      IF (enein(1,1) .EQ. 0.0) enein(1:njene,m) = ene(1:njene)
      DO k =1,nbctim
         DO j=1,njene
            enein(j,k) = enein(j,k)*ene_mult
         ENDDO
      ENDDO
      CALL intrp (-1,-1,renein(1,m),enein(1,m),njene,roa,ene,nj)
      eneb(1) = ene(nj)
      DO j = 1,nj
         IF(ene(j) .LE. 0.0)THEN
             PRINT *,'ERROR electron density intrp  non positive value'
             PRINT *,'at grid point j =',j
             PRINT *,'ene(j) =',ene(j)
             CALL STOP ('subroutine INIT: ene input problem', 0)
         ENDIF
      ENDDO
         CALL chksplnt(nbctim,enein,ksplin,                                     &
                 kbctim,njene,renein,nout,ncrt,ierr,knotsene,'ne')
         CALL prtspln (nbctim,    enein,   renein, bctime, bparene,             &
                     ksplin, kbctim, nout, knotsene, 'ne      ')
         ene_set = .TRUE.
   ELSEIF(inenez .NE.  0 .AND. ene_set == .FALSE.)THEN    ! njene <=0,  ==> make a parabolic profile
          IF ( time .NE. bctime(1))  CALL prblcin (time,bctime,nbctim,enec,     &
                                     eneb,alpene,gamene,ierr,nout,ncrt)
          CALL makpro (r,ene,nj,enec,eneb,alpene,gamene)
          ene_set = .TRUE.
   ENDIF
!
! --- END of initial ne profile input




      IF (njte .GT. 0)THEN                        !spline input
!------------------------------------------------------------------------
! --- te spline profile (at initial time):
!------------------------------------------------------------------------
         ierr1 = 1
         IF (ABS(rtein(2,m)) .LT. pos_def .AND. (jnormin .GE. 3) .AND.  te_set == .FALSE.) THEN 
            !rtein is not set at time m, use rnormin instead:
            CALL copya (rnormin,rtein(1,m),jnormin)
            njte  = jnormin
            ierr1 = 0
         ELSE IF ( ABS(rtein(2,m)) .LT. pos_def .AND. (jnormin .LT. 3) .AND. te_set == .FALSE. ) THEN
            ierr1 = 1
         ELSE
            CALL knotnum(rtein(1,m),njte,ksplin)
            IF (njte .GE. 3)  ierr1 = 0
         END IF
         IF (ierr1 .EQ. 1) THEN
            WRITE  (ncrt, 10) 'te'
            WRITE  (nout, 10) 'te'
10          FORMAT (' ERROR in knots specification of ', a2)
         END IF
         ierr     = MAX0 (ierr, ierr1)
         iprofnbr = kprim + kimp + 2
         IF (ABS(tein(1,1)) .LT. pos_def)  CALL copya (te,tein(1,m),njte)
      DO k =1,nbctim
         DO j=1,njte
            tein(j,k) = tein(j,k)*te_mult
         ENDDO
      ENDDO
         CALL intrp (-1,-1,rtein(1,m),tein(1,m),njte,roa,te,nj)
         DO j = 1,nj
            IF(te(j) .LE. 0.0)THEN
               PRINT *,'ERROR te interpolation yields non positive value'
               PRINT *,'at grid point j =',j
               PRINT *,'te(j) =',te(j)
               CALL STOP ('subroutine INIT: te input problem', 0)
            ENDIF
         ENDDO
         teb(1)   = te(nj)
         CALL chksplnt(nbctim,tein,ksplin,kbctim,njte,rtein,nout,ncrt,      &
                   ierr,knotste,'te')
         CALL prtspln (nbctim,     tein,    rtein, bctime,  bparte,         &
                       ksplin, kbctim, nout,  knotste, 'te      ')
         te_set = .TRUE.
      ELSEIF(te_set == .FALSE.) THEN                                     !parabolic TE input
         IF (time .NE. bctime(1))                                           &
         CALL prblcin (time,bctime,nbctim,tec,teb,                          &
              alpte,gamte,ierr,nout,ncrt)    
         CALL makpro (r,te,nj,tec,teb,alpte,gamte)
         TE_set = .TRUE.
      ENDIF
!
! --- END of initial te profile input
!

! ----------------------------------------------------------------------------
! --- initial  ti profile:
! ----------------------------------------------------------------------------
     IF (njti .GT. 0)THEN
      ierr1 = 1
      IF (rtiin(2,m) .EQ. 0.0 .AND. (jnormin .GE. 3) .AND.  ti_set == .FALSE.) THEN
        CALL copya (rnormin,rtiin(1,m),jnormin)
        njti  = jnormin
        ierr1 = 0
      ELSE IF (rtiin(2,m) .EQ. 0.0 .AND. (jnormin .LT. 3).AND.  ti_set == .FALSE. ) THEN
        ierr1 = 1
      ELSE
        CALL knotnum(rtiin(1,m),njti,ksplin)
        IF (njti .GE. 3)  ierr1 = 0
      END IF
      IF (ierr1 .EQ. 1) THEN
        WRITE (ncrt, 40) 'ti'
        WRITE (nout, 40) 'ti'
      END IF
      ierr     = MAX0 (ierr, ierr1)
      iprofnbr = kprim + kimp + 3
      IF (tiin(1,1) .EQ. 0.0)  CALL copya (ti,tiin(1,m),njti)

      DO k =1,nbctim
         DO j=1,njti
            tiin(j,k) = tiin(j,k)*ti_mult
         ENDDO
      ENDDO
      CALL intrp (-1,-1,rtiin(1,m),tiin(1,m),njti,roa,ti,nj)
      tib(1)   = ti(nj)
      DO j = 1,nj
         IF(ti(j) .LE. 0.0)THEN
             PRINT *,'ERROR ti interpolation yields non positive value'
             PRINT *,'at grid point j =',j
             PRINT *,'ti(j) =',ti(j)
             CALL STOP ('subroutine INIT: ti input problem', 0)
         ENDIF
      ENDDO
      CALL chksplnt(nbctim,tiin,ksplin, kbctim,njti,rtiin,nout,ncrt,         &
                        ierr,knotsti,'ti')
      CALL  prtspln (nbctim,     tiin,    rtiin, bctime,  bparti,            &
                     ksplin, kbctim, nout,  knotsti, 'ti      ')
         ti_set =.TRUE.
      ELSEIF(ti_set == .FALSE.)THEN
         IF (time .NE. bctime(1))  CALL prblcin (time,bctime,nbctim,tic,     &
                                    tib,alpti,gamti,ierr,nout,ncrt)
         CALL makpro (r,ti,nj,tic,tib,alpti,gamti)
         ti_set =.TRUE.

      ENDIF
!
! --- END of initial ti profile


! -------------------------------------------------------------------------------
! --- start initial current spline input:
!--------------------------------------------------------------------------------
    IF (njcur .GT. 0)THEN
      ierr1 = 1
      IF (rcurdein(2,m) .EQ. 0.0 .AND. (jnormin .GE. 3) .AND. curden_set == .FALSE.) THEN
        CALL copya (rnormin,rcurdein(1,m),jnormin)
        njcur = jnormin
        ierr1 = 0
      ELSE IF (rcurdein(2,m) .EQ. 0.0 .AND. (jnormin .LT. 3).AND. curden_set == .FALSE.) THEN
         ierr1 = 1
      ELSE
        CALL knotnum (rcurdein(1,m), njcur, ksplin)
        IF (njcur .GE. 3)  ierr1 = 0
      END IF
      IF (ierr1 .EQ. 1) THEN
        WRITE (ncrt, 40) 'curden'
        WRITE (nout, 40) 'curden'
      END IF
      ierr     = MAX0 (ierr, ierr1)
      iprofnbr = kprim + kimp + 5
      IF (curdenin(1,1) .EQ. 0.0)                              &
          CALL copya (curden, curdenin(1,m), njcur)
      CALL copya (curdenin(1,m),profin,njcur)
      CALL intrp (-1,-1,rcurdein(1,m),profin,njcur,roa, curden, nj)
      CALL chksplnt(nbctim,curdenin,ksplin,kbctim,njcur,rcurdein,nout,ncrt,      &
                                          ierr,knotscur,'curden')
      CALL prtspln (nbctim, curdenin, rcurdein, bctime, bparcur,                 &
                     ksplin, kbctim, nout, knotscur, 'curden  ')
      curden_set = .TRUE.
    ELSEIF(curden_set == .FALSE.)THEN
      IF ( time .NE. bctime(1))  CALL prblcin (time,bctime,nbctim,xjc,           &
                                      xjb,alpxj,gamxj,ierr,nout,ncrt)
      CALL makpro (r,curden,nj,xjc,xjb,alpxj,gamxj)
      curden_set = .TRUE.
    ENDIF
!
! ----------------------------------------------------------------------
! inital angular rotation speed profile (rad/sec)
! ----------------------------------------------------------------------
!
    IF (iangrot .NE. 0)THEN
       IF (ABS(rangrot(1,m))  .GT. pos_def)THEN
          WRITE  (nout, 20)  m, rangrot(1,m)
          WRITE  (ncrt, 20)  m, rangrot(1,m)
20        FORMAT (' ERROR: rangrot(1,',i2,') = ', e14.8, ' must be 0.0')
          ierr = 1
          rangrot(1,m) = zeroc
       ENDIF
       CALL knotnum(rangrot(1,m),njinang,ksplin)

       !
       ! --- rangrot(j,m) is not set correctly, try to copy rnormin
       !
       IF (jnormin .GE. 3 .AND. njinang ==0 .AND. angrot_set == .FALSE. ) THEN
          CALL copya (rnormin, rangrot(1,m), jnormin)
          njinang = jnormin
       ELSEIF(njinang == 0 .AND. angrot_set == .FALSE.)THEN
          WRITE  (nout, 30)  ksplin, m
          WRITE  (ncrt, 30)  ksplin, m
30        FORMAT (' ERROR: array rangrot(1..',i2,',',i2,') is not set'   /     &
                    '        correctly. Last value at any time point must' /   &
                    '        equal 1.0 exactly.')
          ierr    = 1
          njinang = 3
       END IF

       iprofnbr = kprim + kimp + 6
       IF(ABS(angrotin(1,m)) .LT. pos_def)angrotin(1:nj,m)=angrot(1:nj)
      DO k =1,nbctim
         DO j=1,njinang
            angrotin(j,k) = angrotin(j,k)*ang_mult
         ENDDO
      ENDDO
       CALL intrp (-1, -1, rangrot(1,m), angrotin(1,m), njinang, roa,           &
            angrot, nj)
       CALL  prtspln (nbctim, angrotin,  rangrot, bctime, bparang,              &
                     ksplin, kbctim, nout, knotsang, 'omega   ')
       angrot_set = .TRUE.
       !
       ! --- no parabolic profiles for angular rotation:
       !
    ENDIF

    CALL chksplnt(nbctim,angrotin,ksplin,                               &
                 kbctim,njinang,rangrot,nout,ncrt,ierr,knotsang,'omega')

!---------------------------------------------------------------------
! --- initial zeff input:
! --- inenez .NE. 0 ==> zeff from input
!----------------------------------------------------------------------
    IF (njzef .NE. 0 .AND. inenez .NE.  0 )THEN
      ierr1 = 1
      IF (rzeffin(2,m) .EQ. 0.0 .AND. (jnormin .GE. 3).AND. zeff_set == .FALSE.) THEN
        CALL copya (rnormin,rzeffin(1,m),jnormin)
        njzef = jnormin
        ierr1 = 0
      ELSE IF (rzeffin(2,m) .EQ. 0.0 .AND. (jnormin .LT. 3).AND. zeff_set == .FALSE.) THEN
        ierr1 = 1
      ELSE
        CALL knotnum (rzeffin(1,m), njzef, ksplin)
        IF (njzef .GE. 3)  ierr1 = 0
      END IF
      IF (ierr1 .EQ. 1) THEN
        WRITE (ncrt, 40) 'zeff'
        WRITE (nout, 40) 'zeff'
      END IF
      ierr     = MAX0 (ierr, ierr1)
      iprofnbr = kprim + kimp + 4
      IF (zeffin(1,1) .EQ. 0.0)  CALL copya (zeff,zeffin(1,m),njzef)


      DO k =1,nbctim
         DO j=1,njzef
            zeffin(j,k) = zeffin(j,k)*zeff_mult
         ENDDO
      ENDDO
      CALL intrp (-1,-1,rzeffin(1,m),zeffin(1,m),njzef,roa, zeff, nj)
      zeffb(1) = zeff(nj)
      CALL chksplnt(nbctim,zeffin,ksplin, kbctim,njzef,rzeffin,nout,ncrt,     &
                          ierr,knotszef,'zeff')
      CALL prtspln (nbctim,   zeffin, rzeffin, bctime, bparzeff,              &
                          ksplin, kbctim, nout, knotszef, 'zeff    ')
      zeff_set = .TRUE.
    ELSEIF(inenez .NE. 0.AND. zeff_set == .FALSE.)THEN
      IF (time .NE. bctime(1))  CALL prblcin (time,bctime,nbctim,zeffc,       &
                                zeffb,alpzef,gamzef,ierr,nout,ncrt)        
      CALL makpro (r,zeff,nj,zeffc,zeffb,alpzef,gamzef)
      zeff_set = .TRUE.
    ENDIF
!
! --- END initial zeff input
!


 
!---------------------------------------------------------------------------
! --- primary and impurity densities
! --- inenez = 0, calculate electron density and zeff from primary 
! --- and impurity ion densities.
! --------------------------------------------------------------------------
      IF(inenez == 0)THEN
         DO  i=1,nprim
            IF (njenp(i) .NE. 0) THEN            ! input is of spline type
               IF (ABS(renpin(2,i)- 0.0_DP) .LT. pos_def .AND. enp_set == .FALSE. ) THEN
                  DO j=1,njin
                     renpin(j,i) = rnormin(j)
                  END DO
               END IF
               iprofnbr = i
               DO j=1,njenp(i)
                  enp(j,i) = enp(j,i) * density_mult(i)
               ENDDO
               CALL copya (enp(1,i),profin,njenp(i))
               CALL intrp (-1,-1,renpin(1,i),profin,njenp(i),roa,en(1,i),nj)
               enpb(1,i) = en(nj,i)
            ELSE                                   ! input is of parabolic type
               IF (time .NE. bctime(1) .AND. enp_set == .FALSE. )              &
                    CALL prblcin (time,bctime,nbctim,enpc(1,i),enpb(1,i),      &
                    alpenp(1,i),gamenp(1,i),ierr,nout,ncrt)
               CALL makpro (r,en(1,i),nj,enpc(1,i),enpb(1,i),                  &
                    alpenp(1,i),gamenp(1,i))
            ENDIF
            density_mult(i) = 1.0D0
         ENDDO
         enp_set = .TRUE.
         !
         !
         IF (nimp .GT. 0)THEN
            DO  i=1,nimp
               k = nprim + i
               IF (njeni(i) .GT. 0)THEN
                  IF (reniin(2,i) .EQ. 0.0 .AND. eni_set == .FALSE. ) THEN
                     DO j=1,njin
                        reniin(j,i) = rnormin(j)
                     END DO
                  END IF
                  k = nprim+i
                  DO j=1,njeni(i)
                     eni(j,i)   = eni(j,i) * density_mult(k)
                  ENDDO
                  iprofnbr = kprim+i
                  CALL intrp (-1,-1,reniin(1,i),eni(1,i),njeni(i),roa,en(1,k), nj)
                  enib(1,i) = en(nj,k)
               ELSEIF(eni_set == .FALSE.)THEN
                  IF (time .NE. bctime(1))                                 &
                       CALL prblcin (time,bctime,nbctim,enic(1,i),enib(1,i), &
                       alpeni(1,i),gameni(1,i),ierr,nout,ncrt)
                  CALL makpro (r,en(1,k),nj,enic(1,i),enib(1,i),        &
                       alpeni(1,i),gameni(1,i))
               ENDIF
            ENDDO
            eni_set = .TRUE.
         ENDIF

         ! ------------------------------------------------------------------------
         ! --- electron density,zeff for  inenez  =0
         ! --  enbeam,enalp may initially be zero
         ! -------------------------------------------------------------------------
         IF(ene_set == .FALSE. .OR.  zeff_set == .FALSE.)THEN
            DO j=1,nj
               ene (j) = 0.0
               zeff(j) = 0.0
               DO k=1,nion
                  ene(j) = ene(j)+z(j,k)*en(j,k)
                  zeff(j) = zeff(j)+zsq(j,k)*en(j,k)
               END DO
               ene(j)  = ene(j) + enbeam(j) + 2.0*enalp(j)
               zeff(j) = zeff(j) + enbeam(j) + 4.0 * enalp(j)
               zeff(j) = zeff(j)/ene(j)
            ENDDO
            DO j=1,nbctim                 ! need ene at all bctime times. assume constant
               ! since ion input is time independent.
               ! But note that ene is recalculated in zen as
               ! necessary
               enein(1:nj,j)  = ene(1:nj)
               zeffin(1:nj,j) = zeff(j)
               !renein(1:nj,j) = renpin(1:nj,j)             ! Changed Nov 19,'09 HSJ
               !rzeffin(1:nj,j) = renpin(1:nj,j)
               renein(1:nj,j)  = roa(:)
               rzeffin(1:nj,j) = roa(:)
            ENDDO
            knotszef(:)= nj
            knotsene(:) = nj
         ENDIF
      ENDIF          !    inenez =0

! -----------------------------------------------------------------------------------
! -- inenez =-99,-98
! -- in this case ene was set above (under inenez .NE. 0)
! -- we need to set primary ions from ene and impurity
! -----------------------------------------------------------------------------------
      IF(inenez == -99 .OR. inenez == -98 .AND. initialize_from_statefile == .FALSE. )THEN
         IF (nimp .GT. 0)THEN
            DO  i=1,nimp
               k = nprim + i
               IF (njeni(i) .GT. 0)THEN
                  IF (reniin(2,i) .EQ. 0.0) THEN
                     DO j=1,njin
                        reniin(j,i) = rnormin(j)
                     END DO
                  END IF
                  k = nprim+i
                  iprofnbr = kprim+i
                  CALL intrp (-1,-1,reniin(1,i),eni(1,i),njeni(i),roa,en(1,k), nj)
                  enib(1,i) = en(nj,k)
                  eni(1:nj,i) = en(1:nj,k)
               ELSE
                  IF (time .NE. bctime(1))                                 &
                       CALL prblcin (time,bctime,nbctim,enic(1,i),enib(1,i), &
                       alpeni(1,i),gameni(1,i),ierr,nout,ncrt)
                  CALL makpro (r,en(1,k),nj,enic(1,i),enib(1,i),        &
                       alpeni(1,i),gameni(1,i))
               ENDIF
            ENDDO
         ENDIF
         !-------------------------------------------------------------------------
         !now get primary ion densitites:
         !-------------------------------------------------------------------------
          DO i=1,nprim
             renpin(1:nj,i) = renein(1:nj,1)
          ENDDO
          DO j=1,nj 
             en(j,1) = ene(j)
             DO k=nprim+1,nion
                en(j,1) = en(j,1) - z(j,k)*en(j,k)
             END DO
             en(j,1) = en(j,1) - enbeam(j) - 2.0*enalp(j)
             enp(j,1) = en(j,1)
             IF (nprim .EQ. 2) THEN
                en(j,2)  = (1.0 - zfrac)*en(j,1)
                en(j,1)  = zfrac*en(j,1)
                enp(j,1) = en(j,1)
                enp(j,2) = en(j,2)
             END IF
          END DO

          IF (inenez .EQ. -98) THEN

             DO k=1,nimp
         	
         	IF(namei(k).EQ.'he') THEN
         		DO j=1,nj
         		  en(j,nprim+k) = taupin * dtnfus(j)
         		END DO
         	END IF		  
         
             END DO
      
          END IF

          DO j=1,nj
             zeff(j) = 0.0
             DO k=1,nion
                zeff(j) = zeff(j)+zsq(j,k)*en(j,k)
             END DO
             zeff(j) = zeff(j) + enbeam(j) + 4.0 * enalp(j)
             zeff(j) = zeff(j)/ene(j)
          END DO
        
          DO j =1,nbctim
             DO i =1,nion
                   en_bc_inpt(1:kj,i,j)  = en(1:kj,i)
                   ken_bc_inpt(i,j)      = nj            !used in set_gcnmp_vars
                   ren_bc_inpt(1:kj,i,j) = renpin(1:kj,1)
             ENDDO
          ENDDO
       ELSEIF(inenez == -99 .OR. inenez == -98 .AND. initialize_from_statefile == .TRUE. )THEN
            WRITE  (ncrt, 48)
            WRITE  (nout, 48) 
            ierr =1
48          FORMAT (2x,'inenez =-98,-99 not compatible with initialization from statefile')
       ENDIF



! ------------------------------------------------------------------------
! --- transp beam driven current(at initial time):
! -------------------------------------------------------------------------
      IF (external_beam_cur .GT. 0) THEN
          ierr1 = 1
          IF (rcurb_external(2,m) .EQ. 0.0 .AND. (jnormin .GE. 3)) THEN
             CALL copya (rnormin,rcurb_external(1,m),jnormin)
             njcurb_external  = jnormin
             ierr1 = 0
          ELSE IF (rcurb_external(2,m) .EQ. 0.0                          &
                                  .AND. (jnormin .LT. 3)) THEN
             ierr1 = 1
          ELSE
            CALL knotnum(rcurb_external(1,m),njcurb_external,ksplin)
            IF (njcurb_external .GE. 3)  ierr1 = 0
          END IF


          iprofnbr = kprim + kimp + 8
          IF (curbeam_external(1,1) .EQ. 0.0) ierr1 = 1
          ierr     = MAX0 (ierr, ierr1)
          IF (ierr1 .EQ. 1) THEN
            WRITE  (ncrt, 50) 'curbeam_external'
            WRITE  (nout, 50) 'curbeam_external'
50          FORMAT (' ERROR in knots specification of ', a20)
          END IF
          CALL intrp (-1,-1,rcurb_external(1,m),curbeam_external(1,m),      &
                       njcurb_external,roa,curb_external,nj)
          CALL  prtspln (nbctim,    curbeam_external, rcurb_external,       &
                         bctime, bparcurb_external,                         &
                  ksplin, kbctim, nout, knotscurb_external, 'curbtran')
      ENDIF
!
! --- END of initial transp beam driven current  profile input


      IF(diffeq_methd == 3 .OR. inenez .NE. 1) CALL check_en_bc



! --------------------------------------------------------------------------
! --- qine and qini for rfmode ='input'
! --- currently these values are input in inone only 
! --- note that qrfe and qrfi from statefile  may be used as well
! --------------------------------------------------------------------------
      DO j=1,krf
         tstring = ADJUSTL(rfmode(j))
         CALL to_upper_case1(tstring)
         IF(tstring(1:5) .NE. 'INPUT')CYCLE
         IF (njqin(j) .EQ. 0)THEN
            CALL STOP('INIT,njqin must be set',1)
         ENDIF
         njin = njqin(j)
         CALL copya (qine(1,j),profin,njin)
         CALL intrp (0,1,rnormin,profin,njin,roa,qine(1,j),nj)
         CALL copya (qini(1,k),profin,njin)
         CALL intrp (0,1,rnormin,profin,njin,roa,qini(1,j),nj)
      ENDDO








!     bctime is before time0, we must get the profiles at time0 by interpolation:

      IF (time .GT. bctime(1)) THEN
!
        IF (time .LE. bctime(nbctim)) THEN
!
! --- get profiles at initial time:
!
            IF (njene   .NE. 0)  CALL tsplinew (enein, renein, ene,            &
                                     knotsene, kprim+kimp+1,xdum, ydum)
            IF (njte    .NE. 0)  CALL tsplinew (tein, rtein, te,               &
                                knotste, kprim+kimp+2, xdum, ydum)
            IF (njti    .NE. 0)  CALL tsplinew (tiin, rtiin, ti,               &
                                knotsti, kprim+kimp+3, xdum, ydum)
            IF (njzef   .NE. 0)  CALL tsplinew (zeffin, rzeffin, zeff,         &
                                knotszef, kprim+kimp+4,xdum, ydum)
            IF (njcur   .NE. 0)  CALL tsplinew (curdenin, rcurdein,            &
                                curden, knotscur, kprim+kimp+5,xdum, ydum)
            IF (iangrot .NE. 0)  CALL tsplinew (angrotin, rangrot,             &
                                 angrot, knotsang, kprim+kimp+6,xdum, ydum)

          ELSE
            ierr = 1
            WRITE  (nout, 60)  time, bctime(nbctim)
60          FORMAT (' ERROR: time, bctime(nbctim) =', 2(2x, f12.5))
          END IF
        ELSE IF (time .LT. bctime(1)) THEN
          ierr = 1
          WRITE  (nout, 70)  time, bctime(1)
70        FORMAT (' ERROR: time, bctime(1) =', 2(2x, f12.5))
        END IF


      RETURN
40      FORMAT (' ERROR in knots specification of ', a2)

      END SUBROUTINE set_initial_profiles



 
      SUBROUTINE chksplnt (nbctim, profin, ksplin, kbctim, njinpt,         &
                          rprofin, nout, ncrt, ierr1, knots, profile)
!--------------------------------------------------------------------------
! --- SUBROUTINE checks and sets up time-dependent spline profiles
!--------------------------------------------------------------------------
!
      USE nrtype,                            ONLY : I4B,DP
      IMPLICIT  NONE

      INTEGER(I4b) j,nbctim,ksplin,kbctim,njinpt,nout,ncrt,ierr1,          &
           njin,m,njin1

      CHARACTER *(*) profile

      REAL(DP) profin(ksplin,kbctim), rprofin(ksplin,kbctim)
      INTEGER(I4B) knots(kbctim)
!
      njin = njinpt
!
      DO m=1,nbctim
         IF (         m  .GT. 1  )  THEN
            IF (profin(1,m) .NE. 0.0)  go to 230
            !
            ! --- central value at time m is zero. This means time m-1 is to be used
            !
            DO j=1,njin1                    !njin1 is defined for m > 1
               profin(j,m) =  profin(j,m-1)
               rprofin(j,m) =  rprofin(j,m-1)
            END DO
            !
230         IF (rprofin(2,m) .NE. 0.0)  go to 232
            !
            ! --- second knot at time m is zero. this means time m-1 should be used
            !
            DO j=1,njin1
               rprofin(j,m) = rprofin(j,m-1)
            END DO
232         IF (rprofin(1,m) .EQ. 0.0)  go to 235
            WRITE (nout, 10)  m, profile, (rprofin(j,m), j=1,ksplin)
            WRITE (ncrt, 10)  m, profile, (rprofin(j,m), j=1,ksplin)
            ierr1 = 1
235         CALL knotnum (rprofin(1,m), njin, ksplin)
            IF (njin .LT. 3)THEN
               WRITE (nout, 10)  m, profile, (rprofin(j,m), j=1,ksplin)
               WRITE (ncrt, 10)  m, profile, (rprofin(j,m), j=1,ksplin)
               ierr1    = 1
            ENDIF
         ENDIF
         knots(m) = njin
         njin1    = njin
      END DO
!
      RETURN
!
 10 FORMAT (' ERROR in spline knot specification'      /        &
             ' at time point', i5, ' for profile =', a5 /       &
             ' knot values:'                            /       &
               (5(2x, 1pe10.3)                          /))
!
      END SUBROUTINE chksplnt


