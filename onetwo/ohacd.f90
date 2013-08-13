
  
     SUBROUTINE  ohacd(irfc)
!
! ----------------------------------------------------------------------
!          OHMIC HEATING AND CURRENT DRIVE
! --------------------------------------------------------HSJ--06/1504--
! 
       USE param,ONLY  : kj,krf
       USE soln,ONLY   : rbp,curden,curpar_soln,etor,curtype,              &
                         cur_seed,u
       USE geom,ONLY   : fcap,gcap,hcap
       USE mesh,ONLY   : r,ra,roa,dr
       USE numbrs,ONLY : nj,nk,nion
       USE tfact,ONLY  : jhirsh,jneo,wneo
       USE sourc,ONLY  : curboot,eta,etap,curboot_bp,curboot_bt,xjbte,     &
                         xjbne,xjbti,xjbni,xjbnf,curbe,curbi,curbet,curdri,&
                         curb_external,curdbeam,currf,curohm,qohm,         &
                         wohm,qmag,mult_curboot
       USE tcoef,ONLY  : d,dudr,xdit,xdin,xden,dfion,dfast,xdet
       USE mixcom,ONLY : w2mix
       USE ions,ONLY   : z
       USE nub2,ONLY   : enbeam,ibcur
       USE fusion,ONLY : enalp
       USE soln2d,ONLY : rbsaxis
       USE solcon,ONLY : external_beam_cur,q0_max,q0_radius,cur_seed_amp,  &
                         cur_seed_ctr,cur_seed_width,q0_mult,time
       USE rf,ONLY  : extcurrf
       USE extra,ONLY  : q,voltoh,pohm
       USE constnts,ONLY : pi,pisq,twopi
       USE machin,ONLY   : rmajor
       USE tordlrot,ONLY : iangrot
       USE  bd_condtn,ONLY : vloop_obtained,totcur,iter_vloop,vloop_current, &
                             u_vloop_bc,ub
       USE verbose,ONLY : vloopvb
       USE io,ONLY : ncrt
       IMPLICIT NONE
       INTEGER i,j,k,kfar,nclboot_fail,jrbs,irfc,iostat,iodump
       REAL *8 curboot_prev(kj)
       REAL *8                                                             &
              cee,za1,zaz,denfast,dcoef,djbsc,cjrbs,hcapa,extmult,         &
              curqo,const,curohm_nj,tot_curden,curden_nj,rbp_nj,ub_p,      &
              etor_njm1,curdria
!
       LOGICAL   exists, opened
       CHARACTER timestamp*20
       DATA cee     / 2.99792458e10 /
!
!
       etor(:) = 0.0D0
  
!
!
!
!----------------------------------------------------------------------
!         OHMIC HEATING AND CURRENT DRIVE
!-----------------------------------------------------HSJ-6/15/04------
!this section revised (PAP 10/89)
!
!current density -- curden is <J_phi*R_0/R>, the toroidal current
!         inside a flux surface is integral(2 * pi * r*H*dr*curden)
!      units are A/cm**2, which accounts for the strange coefficients
!      E(volt/cm)  = (3.0e2)*E(statvolt/cm)
!      J(A/cm**2)  = (1.0/3.0e9)*J(statamp/cm**2)
!      eta(ohm-cm) = (9.0e11)*eta(sec)
!                     WHERE 3 = 2.99792458 and 9-8.98755179
!
!
      CALL curcalc (rbp, fcap, hcap, gcap,r, curden,curpar_soln, nj)
!
!
!
!next: bootstrap current, which depends on the electron density
!gradient, so some contortions are needed.  This is parallel current.
!    calculate the bootstrap current at mesh centers
!    (RESULT is changed from statamps/cm**2 to amps/cm**2)
!
      kfar = nk - iangrot
!
!For OLD Houlberg bootstrap model (jhirsh = 95 or  96) or
!for NEW NCLASS model             (jhirsh = 99 or 100)
!
 
!
      BOOTSTRAP_TYPE : IF (jhirsh .EQ. 95 .OR. jhirsh .EQ. 96) THEN
         CALL nclboot (jhirsh, curboot,nclboot_fail)
         IF(nclboot_fail .EQ. 0)THEN
               curboot_prev(:) = curboot(:)
         ELSE
               curboot(:) = curboot_prev(:) !if failed use old bootstrap
         ENDIF
        
         DO   j=1,nj  -1
           DO k=1,kfar-1
             d(kfar,k,j) = 0.0
           END DO
         END DO
         curboot(:) = mult_curboot*curboot(:)
!
      ELSE IF (jhirsh .EQ. 99 .OR. jhirsh .EQ. 100) THEN BOOTSTRAP_TYPE
!
         CALL nclass_dr (wneo, jneo, jhirsh, curboot, curboot_bt,        &
                      curboot_bp, eta)
!
      ELSE BOOTSTRAP_TYPE
!
!    
            !Hirshman 88 and Sauter model
         BOOTSTRAP: DO  j=1,nj-1    
            curboot(j) = 0.0
            !
            !calculations are on the half grid
            !ion density, ion and electron temperature terms
            !
            DO   k=1,kfar-1

                curboot(j) = curboot(j)-ra(j)/(cee*eta(j))*d(kfar,k,j)*      &
                                        dudr(k,j) / 2.99792458e9      !amps/cm**2
            ENDDO
               !
               !fast ions (note: assume Z_beam = 1)
               !
               za1 = 0.5 * (z(j,1)+z(j+1,1))
               denfast = (enbeam(j+1)-enbeam(j)+2.0*enalp(j+1)-2.0*enalp(j))    &
                    /(r(j+1)-r(j))
               !
               !IF Hirshman calculation is performed in DIFFUS the electrons
               !associated WITH the fast ions require a different coefficient
               !in the bootstrap calculation
               !
               IF (jhirsh .NE. 0) THEN
                  dcoef = dfast(j) + dfion(j)
               ELSE
                  dcoef = d(kfar,1,j)/za1
               END IF
               curboot(j) = curboot(j)-ra(j)/(cee*eta(j))*dcoef           &
                    * denfast / 2.99792458e9
               !
               ! calculate individual components of bootstrap current IF Hirshman opted for
               !
               IF (jhirsh .EQ. 0)  CYCLE BOOTSTRAP
               k = nion+1
               xjbte(j) = -ra(j)/(cee*eta(j))*xdet(j)*                  &
                    dudr(k,j) / 2.99792458e9
               xjbne(j) = 0.0
               
               DO i=1,nion
                  zaz = 0.5 * (z(j,i)+z(j+1,i))
                  k   = nion + 2
                  xjbti(j,i) = -ra(j)/(cee*eta(j))*xdit(i,j)*            &
                       dudr(k,j) / 2.99792458e9
                  xjbni(j,i) = -ra(j)/(cee*eta(j))*xdin(i,j)*            &
                       dudr(i,j) / 2.99792458e9
                  xjbne(j) = xjbne(j)-ra(j)/(cee*eta(j))*xden(j)*        &
                       zaz*dudr(i,j) / 2.99792458e9
               END DO
               !
               xjbne(j) = xjbne(j)-ra(j)/(cee*eta(j))*dfast(j)          &
                    *denfast / 2.99792458e9
               xjbnf(j) = -ra(j)/(cee*eta(j))*dfion(j)                  &
                    *denfast / 2.99792458e9
               
         ENDDO BOOTSTRAP
 

      END IF BOOTSTRAP_TYPE      ! Option for Bootstrap Model
!
      IF(rbsaxis .GT. 0.0)THEN
         DO j = 1,nj
            djbsc = roa(j)-rbsaxis
            IF(djbsc .GT. 0.0) THEN
                     jrbs = j-1
                     cjrbs = curboot(jrbs)
                     go to 3005
            ENDIF
         ENDDO
 3005    IF(jrbs .GT. 2)THEN
          DO j = 1,jrbs-1
                curboot(j) = cjrbs*((roa(j)/roa(jrbs))**2)
          ENDDO
         ENDIF
      ENDIF
!
!
!     get the driven curren (beam + rf)
!     Note that the beam and rf driven current exist on the full grid
!
      CURRDRIVE : DO  j=1,nj
         IF (ibcur .EQ. 1) THEN
           IF (j .EQ. 1) THEN
               CALL cubicextrp(curbi(2),curbi(3),curbi(4),             &
                    r(2),r(3),r(4),curbi(1),3) 
               CALL cubicextrp(curbe(2),curbe(3),curbe(4),             &
                    r(2),r(3),r(4),curbe(1),3)
               CALL cubicextrp(curbet(2),curbet(3),curbet(4),          &
                    r(2),r(3),r(4),curbet(1),3)
           END IF
           IF(external_beam_cur .EQ. 0)THEN
               curdri  (j) = curbi(j)+curbe(j)+curbet(j)
               curdbeam(j) = curdri(j) ! save total beam driven current
           ELSE
              curdri(j)  = curb_external(j)
              curdbeam(j) =curb_external(j)
           ENDIF
         END IF
!
         extmult = 0.0
!
         DO i=1,krf
           IF (extcurrf(i) .NE. 0.0)  extmult = 1.0
         END DO
!
         IF (irfc .GE. 1 .OR. extmult .NE. 0.0) THEN
           IF (j .EQ. 1) THEN
             CALL cubicextrp (currf(2), currf(3), currf(4),          &
                             r(2), r(3), r(4), currf(1), 3)
           END IF
           curdri(j) = curdri(j) + currf(j)
         END IF
!
      ENDDO CURRDRIVE
!
!
!
!
!
!
!
!
!
!etor: toroidal E field = H*<E.B>/B_T0, the loop voltage is
!   2 * pi * R_0*etor.
!   Ohm's law gives etor = H*eta*(<J.B>/B_T0-J_boot-J_drive)
!         (eta = eta_parallel)
!   units are v/cm.
!   save eta in etap.  eta units are gaussian (sec)
!   neglect toroidal rotation for now
!   at mesh centers 
!
!
      DO j=1,nj-1
        hcapa = 0.5 * (hcap(j)+hcap(j+1))
        curdria = (curdri(j)+curdri(j+1))*0.5
!       Electric field  is determined by the the ohmic current.
!       Note that hcapa*ra(j)*d(kfar,kfar,j)*dudr(kfar,j) is
!       hcapa*eta*<J dot B/Bt0> , d(kfar,kfar,j) = (eta*c^2/(4*pi*F^2*H*r^2)  
!       so the following says 
!       etor = H*eta *(<J .B/Bt0>- <Jboot . B/Bt0> - <Jdrive . B/Bt0>) 
!       or, equivalently,  etor = H*eta * < Johmic . B/Bt0> 
!      here d has units of 1/sec, eta is in sec, dudr is in gauss
!       1.e-8 converts from gauss cm/sec to volts/cm 
        etor(j) = hcapa*(1.0e-8*ra(j)*d(kfar,kfar,j)*dudr(kfar,j)     &
                       -8.98755179e11*eta(j)*(curboot(j)+curdria))
!        curohm(j) = etor(j)/(8.98755179e11*etap(j)*hcapa)
        etap(j) = eta(j)
      END DO
      etor_njm1 = etor(nj-1)
!
! convert curboot, etor, and eta to mesh points
!
      CALL mescon (xjbne,dr,nj)
      CALL mescon (xjbnf,dr,nj)
      CALL mescon (xjbte,dr,nj)
      DO i=1,nion
        CALL mescon (xjbni(1,i),dr,nj)
        CALL mescon (xjbti(1,i),dr,nj)
      END DO
      CALL mescon (curboot,dr,nj)
      CALL mescon (etor,dr,nj)
      CALL mescon (etap,dr,nj)  
!
!--- make sure bootstrap is zero at magnetic axis ... HSJ
!
      curboot(1) = 0.0
!
!
!
!
!
!     ohmic current on  full grid. 
      DO j =1 ,nj
         IF (curtype .EQ. 0) THEN  
!             This is the original definition which doesnt give E = eta*H *johmic
!             because it is based on < Jphi R0/R> instead off < J . B/Bt0>:
              curohm(j)   = curden(j) - curboot(j) - curdri(j)
         ELSE
!             use the parallel definition. Note that this will not give the
!             correct ohmic heating however which is Qohm = etor* < Jphi R0/R>  HSJ
              curohm(j) = etor(j)/(8.98755179e11*etap(j)*hcap(j))
         END IF
      ENDDO
!
!
!
!
!
!
!       prevent q on axis from getting too large, HSJ 08/22/02:
!-------------------------------------------------------------------
         IF(0.0 .LT. q0_max .AND. q0_max .LT. ABS(q(1)))THEN
            !if dq/dt on axis is increasing then etor has
            !negative slope near magnetic axis (on axis slope of etor
            !must be zero however)
            IF((etor(3) - etor(1))/r(3) .LT. 0.0)THEN
               !slope of etor indicates that q is increasing.
               !set slope to zero from r(1) out to r(j)/r(nj)
               cur_seed(:) =0.0d0
            DO j= 2,nj  
                IF(r(j)/r(nj) .LE. q0_radius)THEN
               cur_seed(j) = q0_mult*(etor(j)-etor(1))/(8.98755179e11*etap(j)*hcap(j))
               !add cur_seed to driven current:
               curdri(j) = curdri(j) + cur_seed(j)
               !add cur_seed contribution to etor:
               etor(j)   = etor(j) - 8.98755179e11*etap(j)*hcap(j)*cur_seed(j)
               !redefine curohm to account for cur_seed:
               IF (curtype .EQ. 0) THEN
                  curohm(j)   = curden(j) - curboot(j) - curdri(j)
               ELSE
!                  curohm(j)   = curpar_soln(j) - curboot(j) - curdri(j)
!
!--- use parallel total current rather than toroidal curden to
!--- calculate parallel ohmic current, HSJ:
!
                curohm(j) = etor(j)/(8.98755179e11*etap(j)*hcap(j))
               END IF
              ENDIF
            ENDDO
            ENDIF
            DO j=1,nj,5
              PRINT *,'j, seedcur',j,cur_seed(j)
            ENDDO
!           call exit
         ENDIF
!
 
!
!
!
!
!
!
! calculate the electron ohmic heating term qohm, as E.(J-J(fast ion))
! convert units from W/cm**3 to keV/cm**3-s.
! voltoh: one-turn voltage from ohmic power
!
      DO j=1,nj
        curqo   = curden(j)
        IF (ibcur .EQ. 1) curqo   = curden(j) - curdri(j)
!        qohm(j) = wohm * 0.62415064e16 * etor(j) * curqo
         qohm(j) = wohm * 0.62415064e16 * etor(j) * curohm(j)
        IF (w2mix .LE. 0.0)  qohm(j) = qohm(j) + qmag(j)
      END DO
!
      const  = 4.0 * pisq * rmajor * 1.60217733e-16
      CALL trapv (r,qohm,hcap,nj,pohm)
      pohm   = pohm*const
      voltoh = pohm/totcur(1)
      vloop_obtained = etor(nj)*(twopi*rmajor)        !volts
      
!
!
!
!
!
      tot_curden = totcur(1) !changed below if u_vloop_bc = .TRUE.
      rbp_nj = rbp(nj)
      ub_p = ub(nk-iangrot)
      IF(u_vloop_bc)THEN
         iter_vloop = iter_vloop + 1
         CALL get_vloop_bc(time,vloop_current)
         etor(nj) = vloop_current/(twopi*rmajor)       ! by definition
         vloop_obtained = vloop_current
         !use vloop_current/(2pi*rmajor) = Etor(nj) = eta(nj)*Hcap(nj)*<Johmic B/BT0>(j)
         ! to get curohm(nj) consistent with curret temperatures,etc:
         curohm_nj = curohm(nj)
         curohm(nj) =  vloop_current/(2.*pi*rmajor*8.98755179e11*etap(nj)*hcap(nj))
         !total edge current consistent with vloop_current
         curden_nj = curden(nj)
         curden(nj) = curden(nj) -curohm_nj + curohm(nj)  !amps/cm**2
         !integrate to get  totcur:
         CALL trapv (r, curden , hcap, nj, tot_curden)
         tot_curden = tot_curden*twopi
         rbp(nj) = 0.2*tot_curden
         rbp(nj) = (hcap(nj)*r(nj) + hcap(nj-1)*r(nj-1))*0.5*d(kfar,kfar,nj-1)
         rbp(nj) = dr(nj-1)*0.5*(vloop_current/(twopi*rmajor)    &
                   + etor(nj-1) +0.0*etor_njm1 )/rbp(nj)
         rbp(nj) = -rbp(nj)/300. +  rbp(nj-1)
         rbp(nj) = 0.2*tot_curden
         rbp(nj) =0.5*(rbp(nj)  +rbp_nj)
         u(nk-iangrot,nj) = rbp(nj)
         dudr(nk-iangrot,nj-1) =  (u(nk-iangrot,nj) -  u(nk-iangrot,nj-1))/dr(nj-1)
         ub(nk-iangrot) = rbp(nj)
       ENDIF
!
       !create this file even if u_vloop_bc  = .FALSE.
       !so that we can get the output
         IF( vloopvb .GT. 0)THEN
            iodump = 15
            INQUIRE (file = 'vloop_monitor.txt', iostat = iostat,            &
                                 exist = exists, opened = opened)
            IF (iostat .NE. 0) THEN         ! problem with INQUIRE
               WRITE (ncrt, '(/ 2a, i6)')                                    &
                       ' ERROR: Fatal INQUIRE failure, IOSTAT =', iostat
               CALL STOP ('subroutine ohacd : bad INQUIRE', 12)
            END IF
!
            IF (exists .AND. .NOT. opened) THEN              ! file exists from previous case, trash it
                 CALL DESTROY ('vloop_monitor.txt')
                 exists = .FALSE.
            END IF
            IF(.NOT. exists .AND. .NOT. opened)THEN                             !create the file
                 CALL getioun(iodump,iodump)
                 OPEN (unit = iodump, file = 'vloop_monitor.txt',                          &
                             status = 'NEW', iostat = iostat)
                 IF (iostat .NE. 0) THEN
                    WRITE (ncrt, '(/ 2a, i6)')                                             &
                                 ' ERROR: Fatal OPEN failure, IOSTAT =', iostat
                    CALL giveupus(iodump)
                    CALL STOP ('subroutine DUMP_DATA: bad OPEN', 165)
                 END IF
                 CALL GET_DATE_TIME (timestamp)
                 WRITE (iodump, '(2a /)')                                                  &
                                   'FILE vloop_monitor.txt  CREATED  ', timestamp
                 exists = .TRUE.
                 opened = .TRUE.
            END IF
!
            IF(opened)THEN
                   WRITE (iodump,'(" time = ", 1pe14.6)')time
                   WRITE (iodump,'(" iteration # :",i5," vloop_obtained :",1pe14.6,        &
                     &   " vloop bc :",1pe14.6)')iter_vloop,vloop_obtained,vloop_current
                   WRITE (iodump,'(" total current before correction: ",1pe14.6)') totcur(1)
                   WRITE (iodump,'(" total current after correction: ",1pe14.6)') tot_curden
                   WRITE (iodump,'(" old current density (a/cm**2) at ra  :",1pe14.6 )')curden_nj
                   WRITE (iodump,'(" new current density (a/cm**2) at ra  :",1pe14.6 )')curden(nj)
                   WRITE (iodump,'(" old ohmic current density (a/cm**2) at ra  :",1pe14.6 )')curohm_nj
                   WRITE (iodump,'(" new ohmic current density (a/cm**2) at ra  :",1pe14.6 )')curohm(nj)
                   WRITE (iodump,'(" etor (V/cm) at ra-1 :",1pe14.6)')etor_njm1
                   WRITE (iodump,'(" etor (V/cm) at ra :",1pe14.6)')etor(nj)
                   WRITE (iodump,'(" old  rbp (gauss cm)  at ra :",1pe14.6)')rbp_nj
                   WRITE (iodump,'(" new  rbp (gauss cm)  at ra :",1pe14.6)')rbp(nj)
                   WRITE (iodump,'(" old  ub gauss cm)  at ra :",1pe14.6)')ub_p
                   WRITE (iodump,'(" new  ub gauss cm)  at ra :",1pe14.6)')ub(nk-iangrot)
            ELSE
                   WRITE(ncrt,'(" ERROR in sub ohacd writting file")')
                   CALL STOP('subroutine ohacd :',1)
            ENDIF
         ENDIF
      IF(u_vloop_bc) THEN
          totcur(1) = tot_curden  
      ENDIF
      !totcur(1) will be used to set ub in sub set_boundary_condition
!
      RETURN
!
!
    END SUBROUTINE ohacd
