C******************** START FILE NBIGC.FOR ; GROUP NBIGC ******************
C--------------------------------------------------------------------
C  NBIGC
C  EVALUATE AUXILLIARY PARAMETERS ASSOCIATED WITH THE *GUIDING CENTER*
C  OF A FAST ION ORBIT
C
C  RGA OCT00 - removed levgb = 1,2
C
C--------------------------------------------------------------------
C     PASSED ARGUMENT  INIT=1 IF GC has moved and velocities changed
C                             (reevaluate everything)
C
C                      INIT=0 IF only velocities changed
C                             (reevaluate vperp, vpll, larmor radius)
C
C                      INIT=-1 IF only zone indices are needed.
C
C--------------------------------------------------------------------
C
C  INPUT:  BASIC ORBIT PARAMETERS
C     * THESE QUANTITIES ARE SAVED IN THE MONTE CARLO PTCL LIST ARRAYS
C       BTWEEN BEAM CODE TIMESTEPS
C
C   POSITION:  XION  (REL. FLUX RADIAL COORDINATE)
C                TH  (GENERALIZED POLOIDAL ANGLE)
C   VELOCITY:  VION  ION SPEED (MOD(V))
C             XKSID  ION PITCH (V PARALLEL / VION )
C
C   Also INPUT:
C             BIMODB(MAXION)
C                    MOD(B) AT THE GUIDING CENTER
C                    IF =0 THIS IS RE-EVALUATED HERE
C
C  OUTPUT:  AUXILLIARY GUIDING CENTER ORBIT PARAMETERS
C
C   POSITION:  RMJION  DISTANCE, GUIDING CENTER TO MACHINE AXIS
C                YION  DISTANCE ABOVE/BELOW PLASMA MIDPLANE
C              XBJACO(2,2)  2D JACOBIAN D(R,Y)/(X,TH) *CONDITIONAL*
C
C                RLAR  LARMOR RADIUS OF GYRO
C              CSTHTK  SERIES OF COS(J*TH)  ** GENL CODE ONLY **
C              SNTHTK  SERIES OF SIN(J*TH)  ** GENL CODE ONLY **
C
C   VELOCITY:      EP  VION*VION
C
C                VPLL  V PARALLEL = VION * XKSID
C               VPERP  V PERP =VION * SQRT (1 - XKSID**2)
C
C   INDICES:      NGC  TRANSP ZONE INDEX OF GUIDING CENTER
C                NGC2  BEAM 2D GRID ZONE INDEX OF GUIDING CENTER
C
      SUBROUTINE NBIGC(INIT)
      use map_larmor_mod, only: imap_larmor_badgc  ! flag for additional
                                                   ! guiding center data
C
      use xstraln_calls
      use nbi_com
C
      COMMON/ZTBUG2/ ZPOSTA,ZPOSTP,ZPNBGC,ZPNBFL,ZPGOOS
C
C$	INTEGER INBIGC
C$	SAVE INBIGC
C$	DATA INBIGC/0/
C 1/PI ...
      DATA ZPINV/0.31830985/
C--------------------------------------------------------------
C
C$	CALL CHKCPU(1,INBIGC,ZPNBGC)
C
      imap_larmor_badgc=.true.          ! signal to map_larmor:  NBIGC called.
C
C  EVALUATE INDICES
C
C exit if the particle is outside the coordinate system.
C
      if(init.ne.0) then
C
C  GC moved of indices wanted
C
         IF(XION.GE.XBMBND) then
            NGC=0
            RETURN
         ENDIF
C
         NGCX=IFIX(LCENTR+XION*NZONES)
         NGC=MIN0(LEP1,NGCX)
C
         XINGCX=(XION-XIBLO(NGCX))/(XIBLO(NGCX+1)-XIBLO(NGCX))
         XINGCX=AMAX1(0.0,AMIN1(1.0,XINGCX))
         IF(NGC.LE.LEDGE) THEN
            XINGC=XINGCX
         ELSE
            XINGC=0.0
         ENDIF
C
         IBR=IFIX(LCENTR+NZNBMR*(XION/XMINBM))
         INZ0=NTHZSM(IBR-1)
         INZ=NTHZSM(IBR)
         IF(NLSYM2B) THEN
C  MIDPLANE SYMMETRY ASSUMED; THETA RANGING FROM -PI TO PI;
C  GRID RANGES FROM 0 TO PI
            NGC2=MIN0(INZ,IFIX(INZ0+1+(ABS(TH)*ZPINV)*(INZ-INZ0)))
         ELSE
C  ASYMMETRIC THETA GRID
            NGC2=MIN0(INZ,
     >         IFIX(INZ0+1+((TH-THBDY0)*0.5*ZPINV)*(INZ-INZ0)))
         ENDIF
C  DMC - BUGHUNTING
         ztol=2.0e-7*xpi
         IF((TH.LT.-XPI-ztol).OR.(TH.GT.XPI+ztol)) THEN
            WRITE(NONLIN,9901) NGC2,XION,TH,IBR,INZ,INZ0
 
 9901       FORMAT(' ?? NBIGC - INDEX NGC2 = ',I5,', TH OUT OF RANGE'/
     >         '  XI=',1PE13.6,' TH=',1PE13.6,' IBR,INZ,INZ0=',3(1X,I6))
            print *,'xpi =',xpi
            print *,'TH passed through nbi_com.mod'
            CALL bad_exit
         ENDIF
C
         IF(INIT.lt.0) RETURN
      endif
C
C--------------------------------------------------------------
C
      if(init.gt.0) then

C  (R,Y) LOCATION - FOURIER MOMENTS EXPANSION
C  also EVALUATE JACOBIAN
         IJAC=1
         CALL nbMOMRY(XION,TH,IJAC,RMJION,YION,XBJACO) ! splines
ctest         CALL XMOMRY(XION,TH,IJAC,RMJION,YION,XBJACO)
C
      endif
C
C--------------------------------------------------------------
C
C  VELOCITIES
C
      EP=VION*VION
      VPLL=XKSID*VION
      VPERP=SQRT(1.0-XKSID*XKSID)*VION
C
C--------------------------------------------------------------
C
C  LARMOR RADIUS - CHECK GEOMETRY LEVEL
C
C  NON-CIRCULAR CODE: CHECK FIELD
C   MOD(B) WILL BE KNOWN (FROM THE ORBIT ROUTINE) UNLESS THIS IS
C   THE START OF A NEW ORBIT
      IF(BIMODB(1).EQ.0.0) THEN
C  EVALUATE B, 1/B FROM BIFOURIER SPLINE
C  ALSO EVALUATES DERIVATIVES FOR FIRST STEP IN ORBIT EQN.
         BIXI(1)=XION
         BITHST(1)=TH
         CALL BINTRP
         CALL JINTRP			!RTM 26 AUG.1985
      ENDIF
C
      RLAR=VPERP*ABEAM/(9.578E7*XZBEAMI*BIMODB(1))
C  TOROIDAL FIELD AT G.C.
      ZG=FBZ(NGC,1)+XINGC*(FBZ(NGC+1,1)-FBZ(NGC,1))
      BTORI=BZXR*ZG/RMJION
C
C  ION ANGULAR MOMENTUM (MECHANICAL, GYRO AVG)
      BPHION=RMASSB*VPLL*RMJION*BTILTC(2,NGC2)
C
C  ELECTROSTATIC POTENTIAL ENERGY
      EPHION=XZBEAMI*(PHIPRG(NGC)+XINGC*(PHIPRG(NGC+1)-PHIPRG(NGC)))
C
C--------------------------------------------------------------
C
C  EXIT
C
C$	CALL CHKCPU(2,INBIGC,ZPNBGC)
C
C  SAVE LATEST INSIDE-THE-PLASMA COORDINATES, FOR RADIAL CURRENT
C  CALCULATION (DMC 30 SEPT 1994)
C
      IF(XION.LE.1.0) THEN
         XICUR=XION
         THCUR=TH
         RMJCUR=RMJION
         YCUR=YION
         VPLLCUR=VPLL
         VIONCUR=VION
      ENDIF
C
      RETURN
      END
C******************** END FILE NBIGC.FOR ; GROUP NBIGC ******************
