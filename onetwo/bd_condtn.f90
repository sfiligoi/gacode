MODULE bd_condtn
  USE nrtype,                       ONLY : I4B,DP
  USE param,                        ONLY : kbctim,kk,kprim,kimp,kj,ksplin,kj,kion
  USE common_constants,             ONLY : zero_sca => zeroc
  IMPLICIT NONE
  !        Introduce time  dependent ion densities:
  !        REAL *8, DIMENSION(:,:,:),  POINTER :: en_bc_inpt
  !        REAL *8, DIMENSION(:,:,:)  ,ALLOCATABLE  :: ren_bc_inpt 
  !        note that enein and other kinetic profiles are a 
  !         function of time so we dont redefine them here.
  !        INTEGER, DIMENSION(:,:)  ,ALLOCATABLE  :: ken_bc_inpt 
  !       Onetwo uses namelist input so we cant use above arrays:
  REAL(DP) en_bc_inpt(ksplin,kion,kbctim)
  REAL(DP) ren_bc_inpt(ksplin,kion,kbctim)
  REAL(DP) fix_edge_ni_bc_inpt(kion,kbctim)
  INTEGER(I4B) ken_bc_inpt(kion,kbctim) 
  REAL *8, DIMENSION(:,:) ::                                 &
       ucenter(kbctim,kk), uedge(kbctim,kk),                 &
       ualp(kbctim,kk), ugam(kbctim,kk),                     &
       enpc(kbctim,kprim), enpb(kbctim,kprim),               &
       alpenp(kbctim,kprim), gamenp(kbctim,kprim),           &
       enic(kbctim,kimp), enib(kbctim,kimp),                 &
       alpeni(kbctim,kimp), gameni(kbctim,kimp),             &
       qradin(ksplin,kbctim),enein(ksplin, kbctim),          &  
       tein(ksplin,kbctim), tiin(ksplin,kbctim),             &
       curdenin(ksplin,kbctim),zeffin(ksplin,kbctim),        &
       curbeam_external(ksplin,kbctim),vloop_bc(kbctim),     &
       renpin(ksplin,kprim), reniin(ksplin,kimp),            &
       rtein(ksplin,kbctim), rtiin(ksplin,kbctim),           &
       renein(ksplin,kbctim), rzeffin(ksplin,kbctim),        &
       rcurdein(ksplin,kbctim),bc(kbctim,kk),                &
       rcurb_external(ksplin,kbctim),bctime_zone(kj,kk),     &         
       bparenp (4,kprim ), bpareni(4,kimp  ),                &
       bparte  (4,kbctim), bparti (4,kbctim),                &
       bparang (4,kbctim), bparene(4,kbctim),                &
       bparzeff(4,kbctim), bparcur(4,kbctim),                &
       bparkpol(4,kbctim), bparcurb_external(4,kbctim),      &
       eni(kj,kimp),enp(kj,kprim)

  REAL *8, DIMENSION(:) ::                                   &
       tec(kbctim), teb(kbctim), alpte(kbctim),              &
       gamte(kbctim), tic(kbctim), tib(kbctim),              &
       alpti(kbctim), gamti(kbctim), xjc(kbctim),            &
       xjb(kbctim), alpxj(kbctim), gamxj(kbctim),            &
       qradr(ksplin),dnidt(kj,kion),totcur(kbctim),          &
       enec(kbctim), eneb(kbctim),                           &
       alpene(kbctim), gamene(kbctim), zeffc(kbctim),        &
       zeffb(kbctim), alpzef(kbctim), gamzef(kbctim),        &
       rnormin(kj),dnedt(kj), bctime(kbctim),                &
       ub(kk), fluxb(kk),ub_save(kk), ub_rho_edge(kk),       &
       vloop_bc_time(kbctim),density_mult(kion),dtedt(kj),   &
       dtidt(kj),dpidt(kj,kion),dpedt(kj)
  !
  !
  !
  REAL *8      ::                                            &
       ttweak, fusnin, taupin, ticin, voltin,                &
       qcin, twkfar, w33min, zeflim,zeff_mult,               &
       vloop_thresh,vloop_obtained,vloop_current,            &
       cbound(32),ped_grad,ped_nebar,ped_temp,               &
       ped_nGr,ped_mode,ped_nratio,ene_mult,                 &
       te_mult,ti_mult,ang_mult
  !
  !
  INTEGER *4                                                 &      
       nqrad, njte, njti, njene, itweak,                     &
       njzef, njcurb_external,iprofnbr, mtimeprf,            &
       knoterror,nvloop,iter_vloop
  !
  INTEGER *4, DIMENSION(:) ::                                &
       knotscurb_external(kbctim), knotste(kbctim),       &
       knotsti(kbctim), knotsene(kbctim),                 &
       knotszef(kbctim), knotscur(kbctim),                &
       LBOUND(32),pedestal_models(3)
  !
  !
  CHARACTER(LEN=8)       splninpt
  CHARACTER(LEN=16)       profiles_bcondspl(kprim+kimp+6)
  LOGICAL u_vloop_bc,u_vloop_bc_time,use_pedestal 
  DATA u_vloop_bc,u_vloop_bc_time,use_pedestal               &
       /.FALSE.,.FALSE.,.FALSE./
  DATA lbound /32*0/      !default parms for pedestal
  DATA cbound /32*0.0/
  DATA pedestal_models /3*1/
  !
CONTAINS
  !
  SUBROUTINE setup_vloop(time0,timmax,ncrt,nbctim)
    INTEGER j,ncrt,nbctim
    REAL *8 time0,timmax
    !      output is vloop_bc(1..nvloop)  and vloop_bc_time(1..nvloop)
    !
    ! ---------------------------------------------------------------------
    ! check vloop_bc input.
    DO j=1,kbctim
       IF(vloop_bc(j) .LT. vloop_thresh)u_vloop_bc = .TRUE.
       IF(vloop_bc_time(j) .LT. vloop_thresh) u_vloop_bc_time =.TRUE.
    ENDDO
    !
    IF(u_vloop_bc_time)THEN !vloop_bc given on vloop_bc_time  time base
       !sort the arrays in ascending order:
       CALL piksr2(kbctim,vloop_bc_time,vloop_bc)
       !
       !get number of elements
       nvloop =0
       DO j=1,kbctim
          IF(vloop_bc(j) .LT. vloop_thresh)nvloop=nvloop+1
       ENDDO

       !check for valid time array:
       IF(vloop_bc_time(1) .GT. time0 .OR.                         &
            vloop_bc_time(nvloop) .LT. timmax .OR.               &
            nvloop .EQ. 0)THEN
          WRITE(ncrt,166)nvloop,(vloop_bc_time(j),j=1,nvloop)
166       FORMAT(2x,'ERROR,number of valid time elements',        &
               ' in vloop_bc_time is ',i6,/,                           &
               ' the values are :',/,(2(1pe12.6,2x)))
          WRITE(ncrt,167)time0,timmax
167       FORMAT(2x,'time0,timmax =',2(1pe12.6))
          CALL STOP('vloop_bc_time error',1)
       ENDIF
       !
       IF(bctime(1) .LT. vloop_bc_time(1) .OR.                     &
            bctime(nbctim) .GT. vloop_bc_time(nvloop))THEN
          WRITE(ncrt,168)bctime(1),vloop_bc_time(1),               &
               bctime(nbctim),vloop_bc_time(nvloop)
168       FORMAT(2x,'ERROR, vloop_bc_time must be a '              &
               ' superset of bctime:',/,                               &
               ' bctime(1) =',1pe12.5,/,                               &
               ' vloop_bc_time(1) = ',1pe12.5,/,                       &
               ' bctime(nbctim) = ',1pe12.5,/,                         &
               ' vloop_bc_time(nvloop) = ',1pe12.5)
          CALL STOP('vloop_bc_time error',1)
       ENDIF
       u_vloop_bc = .TRUE.
    ELSE IF(u_vloop_bc)THEN     !vloop_bc given on bctime(1:nbctim)  time base
       !in this case vloop_bc is assumed in correct order correspondig to
       !bctime and at least vloop_bc(1) must be valid:
       IF(vloop_bc(1) .GE. vloop_thresh)THEN
          WRITE(ncrt,165)vloop_bc(1),bctime(1)
165       FORMAT(2x,' Error, vloop_bc(1)= ',1pe12.5,/,              &
               ' bctime(1) = ',1pe12.5 )
          CALL STOP('vloop_bc error',1)
       ELSE
          !vloop_bc(1) is set, if vloop is constant user may not
          !have repeated the values. Adjust for this here:
          nvloop = nbctim
          DO j=2,nbctim
             IF(vloop_bc(j) .GE. vloop_thresh)                      &
                  vloop_bc(j) = vloop_bc(j-1)
          ENDDO
          vloop_bc_time(:) = bctime(:)
       ENDIF
       u_vloop_bc_time =.TRUE.
    ELSE
       vloop_bc(1:kbctim) = zero_sca            ! vloop bc not used
    ENDIF
    !
    !

    RETURN
  END SUBROUTINE setup_vloop
  !
END MODULE  bd_condtn
