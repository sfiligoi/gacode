c-----------------------------------------------------------
c qfm_sub.f
c
c PURPOSE:
c  Manage calls to a local QFM simulation.
c
c INPUT LIST:
c  x_in(1):  q 
c  x_in(2):  s=(r/q)*dq/dr
c  x_in(3):  delta
c  x_in(4):  s_delta=r*d(delta)/dr
c  x_in(5):  kappa
c  x_in(6):  s_kappa=(r/kappa)*d(kappa)/dr
c  x_in(7):  shift=d(R_0)/dr
c  x_in(8):  R_0(r)/a
c  x_in(9):  r/a
c  x_in(10): a*dln(ne)/dr
c  x_in(11): a*dln(ni)/dr
c  x_in(12): a*dln(Te)/dr
c  x_in(13): a*dln(Ti)/dr
c  x_in(14): Ti/Te 
c  x_in(15): betae_unit
c  x_in(16): sqrt(mi/me)
c  x_in(17): (a/cs)*nu_ei
c  x_in(18): (a/cs)*gamma_exb
c  x_in(19): zi2
c  x_in(20): ni2/ne
c  x_in(21): sqrt(mi/mi2)
c  x_in(22): a*dln(ni2)/dr
c  x_in(23): zi3
c  x_in(24): ni3/ne
c  x_in(25): sqrt(mi/mi3)
c  x_in(26): a*dln(ni3)/dr
c  x_in(27): ni/ne [quasineutrality not assumed]
c
c OUTPUT LIST:
c  x_out(1): e/chi_GB
c  x_out(2): chie/chi_GB
c  x_out(3): Di/chi_GB
c  x_out(4): chii/chi_GB
c  x_out(5): Di2/chi_GB
c  x_out(6): chii2/chi_GB
c  x_out(7): Di3/chi_GB
c  x_out(8): chii3/chi_GB
c---------------------------------------------------------
c
      subroutine qfm_sub(x_in,x_out)
c
      implicit none
c
      real*8 x_in(27)
      real*8 x_out(8)
      real*8 eps_r
      real*8 chi_conv
      real*8 sdi
      integer i_print
c
c Original code:
c
      INTEGER FNBI,IN,IMAX,I,J,IW,IX
      REAL*8 CT(0:4),U(5,100)
      COMPLEX*16 ZZ(10),RP(10),W
      PARAMETER (FNBI=20)
      INTEGER DCC,IR,IV,PDEV,IK,IST,ITS,ITL,ITERA,ITC,ISB
      INTEGER SEARCHMODE
      PARAMETER (DCC=38)
      LOGICAL DEF_CORR
      REAL*8 CTE(6,1000),CTI(6,1000),PK(10)
      REAL*8 RAT,RT,ENI,EEI,TE,BTOR,OMEGAC,RHOS,DB,CHINORM,RLT,rno
      REAL*8 A,B,C,DS,ZX,TR,RIW,RISB,chicoeff,chieff,chic2
      REAL*8 EI,TAU,FL,THRD,TVR,STR,XIH
      REAL*8 EN,ENH,PST,RIV,RFL,H,H1,KB,WRAT,WSHIFT
      REAL*8 LBT,TAUI,FTR,EIH,EEH,SH,SH2,FLS,GMMA,FLTOT
      REAL*8 FT,EE,XMIN,YMIN,XMAX,YMAX,DX,NLX
      REAL*8 BQ,EQ,ENQ,Z,GQ,BF,ZE,TAUZ,ZEFF,NQ,NI,G
      REAL*8 CETAIN(32),BETAE,TAUEXP,XX(100)
      REAL*8 E,Q,S,CS,KAPPA,KPPA,RAV,OMA,OMA1,SEARCH
      REAL*8 ALP,ALF,KPC,PR,FTRT,WST,WLN,NG,D1,SI,KIQ,KXQ,WM
      REAL*8 VV(DCC),N,WR,WI,RNEQ,WDE,EPS,WRS,WIS
      REAL*8 SCHI,SCHE,SD,SCHQ,SDQ,EA,HPT(5),GK
      REAL*8 ETE,ETI,ETQ,AZ,AZL,ALA,ALAF,GAV
      REAL*8 VEI,VEF,BTA,COL,COLL,EL,EM,LAMB
      REAL*8 CHI(5),CHE(5),D(5),CHQ(5),DQ(5)
      INTEGER LPRINTIN,NDIM,NEQ,NDISP,IRET
      REAL*8 EIC,EEC,ENC,TAUC,FLC,FTC,EQC,ENQC,BETAEC,TAUZC,QC,SC
      REAL*8 ENHC,LIST,ZFS,KPS,CHIC,R,TEC,diag2,diagz
      REAL*8 WEXB,ROT,zkvot1,zCT,Rnpeak,Rnpeakny,vz,vz1,vz2,vz3
      REAL*8 WZIMAX,TOL,SHPE,SCHEF,DEF,DELM, sumflux,fluxqd,zk,rb2,rb3
      REAL*8 VR(10,10),VI(10,10),fluxm,fluxq,elflux,ionflux,qflux,fluxhi
      REAL*8 LTH,LTE,LN,LNH,LTQ,LNQ,rbou,enb,fluxd,fluxid,anafluxd,fluxh
      COMPLEX*16 HQ,WZ,WZP,zdcoeff,zdcoeff1,zdcoeff2,zdqcoeff1ny
      COMPLEX*16 zdicoeff,zdicoeff1,zdicoeff2,zdiag1,zdqcoeff2ny
      COMPLEX*16 zdqcoeff,zdqcoeff1,zdqcoeff2, W0Z, W0I, WRF
      COMPLEX*16 zhcof1,zhcof2,zhcoeff,zhcof1i,zhcof2i,zhcoeffi 
      COMMON/PARAM/ CT,NG
      COMMON/IRET/ IRET
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XIH,IX,IW
      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/GRAD/ ENC,ENHC,EIC,EEC
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BF,GQ,SI,ZE,TAUZC,ENQC,EQC,KIQ,KXQ,ZFS
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUC,AZL,
     &FTC,FLC,LPRINTIN,NDIM
      COMMON/W/ WM,WR,WI
      COMMON/BETAE/ BETAEC,QC,SC,CS,EM
      COMMON/TEST/ ALP,ALF,KPC,KPS
      COMMON/HQ/ H1,HQ,GAV,WZ,WZP,WZIMAX
      COMMON/SHAFRS/ ALAF
      COMMON/KAPPA/ KPPA,RAV
      COMMON/SCALE/ NLX
      COMMON/IK/ IK,IST,ITC,ITL,ITS,ITERA,TOL
      COMMON/TP/ EA,HPT
      COMMON/WROT/ WEXB,ROT
      COMMON/GRKVOT/ GK
      COMMON/ZV/ VR,VI
      COMMON/ZZ/ ZZ
      COMMON/NEQ/ NEQ
      COMMON/EM/ SHPE,CHIC,SCHEF,DEF
      COMMON/LT/ LTH,LTE,LN,LTQ,LNQ
      COMMON/ISB/ ISB
      COMMON/SEARCHMODE/ SEARCHMODE
      COMMON/OMA/ OMA
      COMMON/WDE/WDE
      COMMON/DELM/ DELM,OMA1
C
      EQUIVALENCE (VV(1),EN ),(VV(2),EI ),(VV(3),EE )
     &           ,(VV(4),TAU ),(VV(5),FL ),(VV(6),FT )
     &           ,(VV(7),PDEV ),(VV(8),PST ),(VV(9),XMIN )
     &           ,(VV(10),XMAX),(VV(11),YMIN),(VV(12),YMAX)
     &           ,(VV(13),DX),(VV(14),RIV),(VV(15),BQ)
     &           ,(VV(16),EQ),(VV(17),ENQ),(VV(18),Z)
     &           ,(VV(19),BETAE),(VV(20),AZ),(VV(21),COL)
     &           ,(VV(22),EL),(VV(23),TE),(VV(24),TAUZ)
     &           ,(VV(25),PR),(VV(26),Q),(VV(27),S)
     &           ,(VV(28),E),(VV(29),N),(VV(30),LIST)
     &           ,(VV(31),RNEQ),(VV(32),ALA),(VV(33),KAPPA)
     &           ,(VV(34),TR),(VV(35),RIW),(VV(36),RISB)
     &           ,(VV(37),BTOR),(VV(38),SEARCH)
C

      integer ieee_handler,itest
      CHARACTER*5 CC(DCC),CC1          !!Names may be 5 characters long
      DATA CC( 1)/'EN'/,CC( 2)/'EI' /,CC( 3)/'EE' /
     &    ,CC( 4)/'TAU' /,CC( 5)/'FL' /,CC( 6)/'FT' /
     &    ,CC( 7)/'PDEV' /,CC( 8)/'PST' /,CC( 9)/'XMIN' /
     &    ,CC(10)/'XMAX' /,CC(11)/'YMIN' /
     &    ,CC(12)/'YMAX' /,CC(13)/'DX' /,CC(14)/'RIV' /
     &    ,CC(15)/'BQ' /,CC(16)/'EQ' /CC(17)/'ENQ' /
     &    ,CC(18)/'Z' /,CC(19)/'BETAE' /,CC(20)/'AZ' /
     &    ,CC(21)/'COL' /,CC(22)/'EL' /,CC(23)/'TE' /
     &    ,CC(24)/'TAUZ' /,CC(25)/'PR' /,CC(26)/'Q' /
     &    ,CC(27)/'S' /,CC(28)/'E' /,CC(29)/'N' /
     &    ,CC(30)/'LIST' /,CC(31)/'RNEQ' /,CC(32)/'ALA' /
     &    ,CC(33)/'KAPPA' /,CC(34)/'TR' /,CC(35)/'RIW' /
     &    ,CC(36)/'RISB' /,CC(37)/'BTOR' /,CC(38)/'SRCH' /
C--------

c      itest=ieee_handler('set','common',2)
c      OPEN(UNIT=FNBI,FILE='compf9.start'
c     &,FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='OLD'
c     &,ERR=09950)
c      GOTO 09951
c09950 STOP 'CANNOT OPEN START FILE'
c09951 CONTINUE
c09914 READ(FNBI,'(A)',END=09912)CC1
c      IF(CC1(1:1).EQ.';')GOTO 09914 !!Check for comment line
c      DO 09913 J=1,DCC
c      IF(CC1.EQ.CC(J))THEN
c      READ(FNBI,*)VV(J)
c      GOTO 09914
c      END IF
c09913 CONTINUE
c      STOP 'INPUT ERROR'
c09912 CONTINUE
c      CLOSE(UNIT=FNBI)
C------

      i_print = 0

      if (i_print .eq. 1) then
00010 DO 00011 J=1,DCC
      WRITE(*,*)CC(J),VV(J)
00011 CONTINUE
      endif

C00014 WRITE(*,*)'New value ? OK means continue'
C      READ(*,'(A)')CC1
C      IF(CC1.EQ.'OK')GOTO 00012
C      DO 00013 J=1,DCC
C      IF(CC1.EQ.CC(J))THEN
C      READ(*,*)VV(J)
C      GOTO 00014
C      END IF
C00013 CONTINUE
C      WRITE(*,*)'Input error'
C      GOTO 00010
00012 CONTINUE
C

c--------------------------------------------------------------
c      INPUT DATA

      eps_r = x_in(9)/x_in(8)      

c     Code will fail if r too small, so for now
c     we need this crude test to prevent failure.
      if (eps_r < 0.05) then
         x_out(:) = 0.0
         return
      endif

c     EN = 2Lne/R
      VV(1)=2d0/x_in(10)/x_in(8) 

c     EI = Lni/LTi
      VV(2)=x_in(13)/x_in(11)

c     EE = Lne/LTe (not used) 
      VV(3)=x_in(12)/x_in(10)  

c     TAU = Te/Ti
      VV(4)=1d0/x_in(14) 

c     FLR parameter (fixed)
      VV(5)=0.04d0     

c     Fraction of trapped electrons
      VV(6)=sqrt(2d0*eps_r/(1d0+eps_r))
 
      if (x_in(20) .gt. 0d0) then
c     n_z/n_e
         VV(15)=x_in(20)
c     EQ = Lnz/LTz (=Lnz/LTi)
         VV(16)=x_in(13)/x_in(22)
c     ENQ = 2 LNz/R
         VV(17) = 2d0/x_in(22)/x_in(8)
c     z
         VV(18)=x_in(19)
c     Az (impurity mass; assume main ions are D)
         VV(20)=2d0/x_in(21)**2
      else
         VV(15)=1d-3
         VV(16)=10.0d0
         VV(17)=10.0d0
         VV(18)=6d0
         VV(20)=12d0
      endif


c     betae
      VV(19)=x_in(15)

c     TAUZ = Te/Tz (assume Ti=Tz)
      VV(24)=1d0/x_in(14) 

c     q
      VV(26)=x_in(1)    

c     s
      VV(27)=x_in(2)   

c     a/R
      VV(28)=1d0/x_in(8)

c------------------------------------------------------------------
c     Extra parameters

      VV(7)=9d0
      VV(8)=1d0
      VV(9)=0.04d0    ! Xmin for loop
      VV(10)=0.15d0   ! Xmax for loop
      VV(11)=-15.0d0
      VV(12)=15.0d0
      VV(13)=0.1d0    ! dx for loop
      VV(14)=5d0      ! FLR-scaling,one point

      VV(21)=0d0
      VV(22)=0d0

c     Te
      VV(23)=1d0

c     a
      VV(25)=1d0

c     ne
      VV(29)=1d0      

      VV(30)=1d0
      VV(31)=9d0     ! number of eqs
      VV(32)=1d0
      VV(33)=1d0
      VV(34)=1d0
      VV(35)=1d0
      VV(36)=2d0     !=1 gives strong ballooning

c     Btor
      VV(37)=1d0     

      VV(38)=2d0
  
c--------------------------------------------------------------

      if (i_print .eq. 1) then
      OPEN(UNIT=1,FILE='DATA1'
     &,FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='UNKNOWN'
     &,ERR=00950)
      OPEN(UNIT=2,FILE='DATA2'
     &,FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='UNKNOWN'
     &,ERR=00950)
      OPEN(UNIT=3,FILE='DATA3'
     &,FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='UNKNOWN'
     &,ERR=00950)
      OPEN(UNIT=4,FILE='DATA4'
     &,FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='UNKNOWN'
     &,ERR=00950)
      endif

      GO TO 00951
c
00950 STOP 'CANNOT OPEN DATA FILE'
c
00951 CONTINUE
c
      IV=INT(RIV)

      if (i_print .eq. 1) WRITE(1,00017) XMIN,XMAX,YMAX,IV

00017 FORMAT(2X,3G11.3,I5)
c
      DO 00099 I=1,6
      DO 00099 J=1,200
      CTE(I,J)=0.
00099 CTI(I,J)=0.
c
      THRD=1.D0/3.D0
      TVR=2.D0*THRD
      FTR=5.D0/3.D0
      STR=7.D0/3.D0
      BTA=1.5
      EPS=0.178
      GK=1.
      ZFS=0.D0
      AZL=AZ
      EM=EL
      ZE=Z
      LBT=EI/EN
      ENI=EN
      EEI=EE
      BF=BQ
      G=1.-Z*BQ
      NI=G*N
      ENH=G*EN/(1.-Z*BQ*EN/ENQ)
      EIC=EI
      EEC=EE
      TAUC=TAU
      TEC=TE
      TAUEXP=TAU
      FLC=DSQRT(FL)
      FTC=FT
      EQC=EQ
      ENQC=ENQ
      ENHC=ENH
      BETAEC=BETAE
      TAUZC=TAUZ
      QC=Q
      SC=S
      NEQ=INT(RNEQ)
      ISB=INT(RISB)
      SEARCHMODE=INT(SEARCH)
c      ALAF=ALA
      NDISP=NEQ
      NDIM=5
      KPPA=KAPPA
      IST=1
      IK=1
      WEXB=0.
      ROT=0.
      EA=E
      IW=INT(RIW)
      R=PR/E
      D1=6.6/(R*BTOR**2)
      OMEGAC=0.957*BTOR*1.D8
      if (i_print .eq. 1) WRITE(*,10014) R,BTOR,E,D1
10014 FORMAT(' R=',G11.3,' BTOR=',G11.3,' E=',G11.3,' D1=',G11.3)
      IX=1
      KB=0.4092/BTOR**2
C
      DEF_CORR =.FALSE.
      IF(FL.LT.0.) DEF_CORR =.TRUE.
      FL=DABS(FL)
c
      DO 00015 I=1,31
00015 CETAIN(I)=0.
      CETAIN(32)=1.D-15
      LPRINTIN=LIST
C
c  -----  VECTOR  PK  IS  DEFINED  ------
c
      PK(1)=TAU
      PK(2)=EI
      PK(3)=EN
      PK(4)=FT
      PK(5)=FL
      PK(6)=BETAE
      PK(7)=KAPPA
c      PK(8)=COL
      PK(8)=EE
      PK(9)=Q
      PK(10)=S
C
      IV=INT(RIV)
      FTRT=FTR/TAU
      TAUI=1./TAU
      RFL=SQRT(FL)
      RAT=EI/EN
      GQ=1.-Z*BQ
      NQ=BQ*N
      NI=GQ*N
      ZEFF=(NI+Z*Z*NQ)/N
C-----
      IMAX=INT((XMAX-XMIN)/DX)
c      IMAX=1
C
      if (i_print .eq.1) write(*,*)'IMAX=',IMAX 
c ----------------------------------------- Main loop -------
      DO 00150 IN=1,IMAX
      PK(IV)=XMIN-DX+IN*DX
      XX(IN)=PK(IV)
      if (i_print .eq. 1) WRITE(*,999)
  999 FORMAT(//)
      TAU=PK(1)
      EI=PK(2)
c      EE=EI
      EN=PK(3)
      FT=PK(4)
      FL=PK(5)
      BETAE=PK(6)
      KAPPA=PK(7)
c      COL=PK(8)
      EE=PK(8)
c     !!!OBS - R/LT fix!!!
c       EE=0.5*EN*7.d0
c      EI=0.5*EN*7. ! Se nedan EI=...!!
c       EQ=0.5*ENQ*7.d0
c       EQ=EI
      EQC=EQ
      Q=PK(9)
      S=PK(10)
C
      TAUC=TAU
c      IF(IV.EQ.1) TAUZC=TAUC
      EIC=EI
      EEC=EE
      ENC=EN
      ENQC=ENQ
c      ENHC=EN
      FTC=FT
      QC=Q
      SC=S
c      TE=TEC*TAU/TAUEXP    !!!!  TI kept constant
      RFL=DSQRT(FL)
c
c      GMMA=BSI0E(FL)
c      FLTOT=1.-GMMA
      IF(.NOT.DEF_CORR) GO TO 01121
C  ***********************************
C  New part with varying correlation length
c      SH2=2.*S-1.+(KAPPA*(S-1.))**2
c      SH=DSQRT(SH2)
c      FLS=(0.7+2.4/(7.14*Q*SH+0.1))*FL
c      RFL=DSQRT(2.D0*FLS/(1.D0+1.D0/TAU))
C  *****************************************
01121 CONTINUE
      FLC=RFL
      BETAEC=BETAE
      KPPA=KAPPA
      COLL=COL
      LN=0.5*R*EN/PR
      LTE=LN/EE
      BF=BQ
      G=1.-Z*BQ
      GQ=1.-Z*BQ
      ENH=GQ*EN/(1.-Z*BQ*EN/ENQ)
      ENHC=ENH
c      EI=0.5*ENH*0.d0
      EIC=EI
      LNH=0.5*R*ENH/PR
      LNQ=0.5*R*ENQ/PR
      LTH=LNH/EI
      LTQ=LNQ/EQ
C
c      ENH=G*EN/(1.-Z*BQ*EN/ENQ)
      IF(ABS(ENH).GE.0.001) GOTO 79
      IF(ENH.LT.0.) ENH=-0.001
      IF(ENH.GE.0.) ENH=0.001
   79 CONTINUE
C
      TAUI=1./TAU
      FTRT=FTR/TAU
      RT=EN/ENI
C      EI=RAT*EN
C      EE=RT*RT*EEI
      EIH=EI-7./3.+FTR*EN
      EEH=EE-7./3.+FTR*EN
      GQ=1.-Z*BQ
C     IF(IV.EQ.2) EE=EI
      ETI=ENH/EI
      ETE=EN/EE
      ETQ=ENQ/EQ
C
      KIQ=ENQ/ENH
      KXQ=EN/ENQ
C
c      BETAE=0.01*KB*N*TE
      NI=GQ*N
      CS=0.311*DSQRT(TE)*10**6
      if (i_print .eq.1) WRITE(*,00126) EN,EI,EE,FL,TAU
00126 FORMAT(2X,'EN=',F8.3,' EI=',F8.3,' EE=',F8.3,' FL=',F8.3,
     &' TAU=',G12.4)
      if (i_print .eq.1) WRITE(*,00299) ENQ,ENH
00299 FORMAT(2X,'ENQ=',G11.3,' ENH=',G11.3)
C     WRITE(*,00127) ENI,EEI,RAT,RT
00127 FORMAT(2X,'ENI=',G10.4,' EEI=',G10.4,' RAT=',G10.4
     1,' RT=',G10.4)
      if (i_print .eq.1) WRITE(*,00129) BETAE
00129 FORMAT(2X,'BETAE=',G11.3)
      if (i_print .eq.1) WRITE(*,00131) Q,S,CS
00131 FORMAT(2X,'q=',G12.4,' S=',G12.4,' CS=',G12.4)
      if (i_print .eq.1) WRITE(*,00132) FT
00132 FORMAT(2X,'FT=',G12.4)
c      WRITE(*,00134) ALA
00134 FORMAT(2X,'ALA=',G11.3)
      if (i_print .eq.1) WRITE(*,00135) FL,KAPPA,S,SH,FLS,RFL
00135 FORMAT('  FL=',G11.3,' KAPPA=',G11.3,' S=',G11.3,' SH=',
     &G11.3,' FLS=',G11.3,' RFL=',G11.3)
C
      WST=DSQRT(FL)*CS/(PR*ABS(LN))
      WDE=ABS(EN)*WST
      WLN=CS/(DSQRT(TAU)*PR*ABS(LN))
      LAMB=15.95-DLOG(DSQRT(N)/TE)
      VEI=(0.09*NI*LAMB/TE**1.5)*10**4
      VEF=VEI/(EPS*WDE)
      VEF=COL*VEF
C
      if (i_print .eq.1) WRITE(*,00130) WST,WLN,VEF,ZEFF,COL
00130 FORMAT(' WST=',G11.3,' WLN=',G11.3,' VEF=',G11.3,' ZEFF=',G11.3,
     &' COL=',G11.3)
C
      H=0.5*ABS(S)/q
      if (i_print .eq.1) WRITE(*,00133) H
00133 FORMAT(2X,'H=',G12.4)
C   -----------------------------------------------
C
C   Calculation of simple local ITG growthrate
C-------------------------------------------------
      A=1.-EN*(1.+10./(3.*TAU))-FL*(1+EI+5.*EN/3.)/TAU
      B=EI-7./3.+5.*EN*(1.+1./TAU)/3.+5.*FL*(1.+EI)/(3.*TAU)
      B=B*EN/TAU
      C=A/(2*(1+FL))
      DS=C*C-B/(1+FL)
      IF(DS.LT.0.) GOTO 140
      WR=C+SQRT(DS)
      WI=0.
      GO TO 160
  140 WR=C
      WI=SQRT(-DS)
      if (i_print .eq.1) WRITE(*,170) WR,WI
  160 CONTINUE
  170 FORMAT(2X,'Local ITG eigenvalue WR=',F7.3,' WI=',F7.3)
C-----------------------------------------------------------------------
C Iteration control variables
      ITC=1
      ITL=100
      TOL=0.05
      IST=1
c----------
C Call to linear solver
c
      CALL disp9t(NDISP,ZZ)
c--------
      if (i_print .eq.1) WRITE(*,174) ISB
  174 FORMAT(' ISB=',I5)
      if (i_print .eq.1) WRITE(*,175) ITC,ITS,ITERA
  175 FORMAT('  ITC=',I5,' ITS=',I5,' ITER=',I5)
      if (i_print .eq.1) WRITE(*,177) WZ,WZP
  177 FORMAT('  WZ=',2G11.3,' WZP=',2G11.3)
      OMA1=OMA
C
      if (i_print .eq.1) WRITE(*,00219) OMA
00219 FORMAT(/,' OMA= ',G11.3,/)
      IR=0
C
C      chicoeff=(EI-2.d0/3.d0-(1.-FT)*10.*EN/(9.*TAU))/EI
C      chieff=0.d0
      DO 00199 I=1,NDISP
      ZX=DIMAG(ZZ(I))
      IF(ZX.LE. 0.001) GOTO 00199
      IR=IR+1
      RP(IR)=ZZ(I)
00199 CONTINUE
C
      if (i_print .eq.1) WRITE(*,310) IR
  310 FORMAT(2X,' IR=',I4)
C
      if (i_print .eq.1) WRITE(*,00134) ALAF
      fluxd=0.d0
      fluxh=0.d0
      fluxhi=0.d0
      fluxq=0.d0
      fluxm=0.d0
      fluxid=0.d0
      fluxqd=0.d0
      anafluxd=0.d0
      diagz=0.d0
      vz=0.d0
      rbou=1.d0
c      rbou=0.25+s*2.d0/3.d0
      enb=EN*rbou
c      rb2=2./3.+5.*s/9.
      rb2=1.D0
c      rb3=1.+0.79*s**2
      rb3=1.D0 
      zk=-ENQ*Z*TAUZ
      W0Z=0.d0
      W0I=0.D0
C  Printout of eigenvalues with different normalisations
      DO 0200 I=1,IR
      W=RP(I)
      WR=DREAL(W)
      WI=DIMAG(W)
      if (i_print .eq. 1) then
         WRITE(*,311) WR*(2.*dsqrt(FL)),WI*(2.*dsqrt(FL)),I
      endif
      WRF=W
c      IF(WR.LE.0.) GOTO 0200
      WRS=EN*DREAL(W)
      WIS=EN*DIMAG(W)
      if (i_print .eq. 1) WRITE(*,411) WRS,WIS,I
      if (i_print .eq. 1) WRITE(1,*)XX(IN),WR,WI
      zdcoeff1=W*EN*(1.d0-enb)-(7.d0/3.d0-EE-5.*enb/3.)*enb
      zdcoeff2=W**2*EN**2-10.*W*EN*enb/3.+5.*enb**2/3.
      zdicoeff1=W*EN*(1.d0*LN/LNH-EN)
     &+((7.d0/3.d0)*LN/LNH-EI-5.*EN/3.)*EN/TAU
     &-FL*(W*EN+(1.+EI)*(LN/LNH)/TAU)*(W*EN+5.*EN/(3.*TAU))
      zdicoeff2=W**2*EN**2+10.*W*EN*EN/(3.*TAU)+5.*EN**2/(3.*TAU**2)
      zdqcoeff1=Z*TAUZ*(ENQ*rb2-1.)*(WRF*zk-5.*ENQ*rb2/3.)
     &-zk*rb2*(2./3.-EQ)
      zdqcoeff1ny=(W+5./(3.*TAU*Z))/((W-2./(TAU*Z))*8.*Q**2*FL)
      zdqcoeff2=WRF**2*zk**2-10.*WRF*zk*ENQ*rb2/3.+5.*ENQ**2*rb2**2./3.
c     &-WRF*zk  
      zdqcoeff2ny=W**2+10.*W*rb2/(3.*TAU*Z)+5.*rb2**2/(3.*TAU**2*Z**2)
      zdiag1=W**2+W*ENQ*rb2*10./(3.*TAU*Z)+5.*(ENQ*rb2)**2/
     &(3.*TAU**2*Z**2)
      diag2=DIMAG(W)*((ABS(W))**2+14.*DREAL(W)*rb2/(3.*TAU*Z)
     &+55.*rb2**2/(9.*TAU**2*Z**2))/EN
      vz1=((2.*EQ/ENQ)*2.*DIMAG(W)*rb2/(TAU*Z*EN))*(DREAL(W)
     &+5.*rb2/(3.*TAU*Z))
      vz2=2.*DIMAG(W)*rb2/EN*((ABS(W))**2+DREAL(W)*10.*rb2/(3.*TAU*Z)
     &+35.*rb2**2/(9.*TAU**2*Z**2))
      vz3=DIMAG(W)*(2.*DREAL(W)*(ABS(W)**2-DREAL(W)*rb2/(3.*TAU*Z)
     &-10.*rb2**2/(3.*TAU**2*Z**2))
     &+(7.*DREAL(W)**2-DIMAG(W)**2/3.+100.*DREAL(W)*rb2/(9.*TAU*Z)
     &-5*rb2**2/(Z**2*TAU**2))/(TAU*Z))
     &/(4.*EN*Q**2*FL*(ABS(W-2.*rb2/(TAU*Z)))**2)
      zdcoeff=FT*zdcoeff1/zdcoeff2
      zdicoeff=zdicoeff1/zdicoeff2
      zdqcoeff=zdqcoeff1/zdqcoeff2+zdqcoeff1ny/zdqcoeff2ny
      zkvot1=(DREAL(WRF)**2+DIMAG(WRF)**2)
     &+14.*DREAL(WRF)/(3.*TAU*Z) +55./(9.*TAU**2*Z**2)
      zCT=2.*(-DREAL(WRF)-5./(3.*TAU*Z))/(Z*zkvot1)
      Rnpeak=-zCT*2.*EQ/ENQ+2.*(1.+2.*zCT/3.)
      fluxd=fluxd-DIMAG(zdcoeff)*WIS**2/(EN*DSQRT(FL))
      fluxid=fluxid-DIMAG(zdicoeff)*WIS**2/(EN*LN*DSQRT(FL)/LNH)
      fluxqd=fluxqd-DIMAG(zdqcoeff)*WIS**2/(EN*LN*DSQRT(FL)/LNQ)
      diagz=diagz+(diag2/(ABS(zdqcoeff2ny)**2))*WIS**2/(EN*DSQRT(FL))
      vz=vz-((vz1+vz2-vz3)/((ABS(zdqcoeff2ny))**2))*WIS**2/
     &(EN*DSQRT(FL))
c     or, divide with zdiag1 above?
c      Rnpeakny=-vz/diagz
      anafluxd=anafluxd+(1.-EN)*WIS**3/(EN*DSQRT(FL)*(WIS**2+WRS**2))
      zhcof1=W*EN/(W*EN-5.*enb/3.)
      zhcof2=(2./3.-EE)/(W*EN)
      zhcoeff=2.*zhcof1*zdcoeff/(3.*FT)-zhcof1*zhcof2
      zhcof1i=W*EN/(W*EN+5.*EN*rb2/(3.*TAU))
      zhcof2i=(2./3.-EI)/(W*EN)
      zhcoeffi=2.*zhcof1i*(zdcoeff+1-FT)/(3.)-zhcof1i*zhcof2i
      fluxh=fluxh-FT*DIMAG(zhcoeff)*WIS**2/(EN*EE*DSQRT(FL)) 
      fluxhi=fluxhi-DIMAG(zhcoeffi)*WIS**2/(EN*EI*DSQRT(FL))
      fluxm=fluxm+(1.d0/EN)*WIS**3/(dsqrt(FL)*(WRS**2+WIS**2))
C      chic2=(WI**3/(WI**2+(WR+5./(3.*TAU))**2))/dsqrt(FL)
C      IF(WI.LE.0.00001) GOTO 0200
C      chieff=chieff+chicoeff*chic2
C      WRAT=DABS(WI/WR)
C      WSHIFT=WR+FTRT
C      WRITE(*,01177) WSHIFT,WRAT
 0200 CONTINUE
c------------------------------------
      Rnpeakny=-vz/diagz
      if (i_print .eq.1) WRITE(*,0202) WDE,WIS,WLN,DIMAG(W),FL
 0202 FORMAT(' WDE=',G11.3,' WIS=',G11.3,' WLN=',G11.3,' WI=',G11.3,
     &' FL=',G11.3)
  311 FORMAT(//,'W norm by WDE(or cs/R)WR=',G11.5,'WI=',G11.5,'I=',I5)
  411 FORMAT(/,2X,'WRS=',G11.3,' WIS=',G11.3,' I=',I5)
01177 FORMAT('WR-FTRT=',G11.3,' WI/WR=',G11.3)
C
C      WRITE(*,*)'chieff=', chieff
      WZ=EN*WZ
      HQ=EN*HQ
      if (i_print .eq.1) WRITE(*,00128) ALF,ALP,WZ,KAPPA
00128 FORMAT(/,2X,' ALF=',G11.3,' ALP=',G11.3,' WZ=',2G11.3,
     &/,' KAPPA=',G11.3)
      if (i_print .eq.1) WRITE(*,312) H1,HQ,GAV,RAV
  312 FORMAT(//,2X,'H1=',G12.4,' HQ=',2G12.4,' GAV=',G12.4,
     &/,' RAV=',G12.4)
      U(1,IX)=TAUI*TE
      U(2,IX)=TE
      U(3,IX)=N
      U(4,IX)=TE/TAUZ
      U(5,IX)=BQ*N
C
      fluxq=ENQ*(fluxd/EN-(1.-Z*BQ)*fluxid/ENH)/(Z*BQ)
      elflux=fluxd/LN
      ionflux=fluxid/LNH*(1.-Z*BQ)
c      qflux=fluxqd/LNQ*(Z*BQ)
      qflux=fluxqd*2./ENQ
      sumflux=elflux-ionflux-qflux
c  Calculation of transport coefficients

c
      CALL DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DQ)
C
c---------------------------------------------
      rno=D1*(Te**1.5)
 319  if (i_print .eq.1) then
         WRITE(*,330) SCHI/rno,SCHE/rno,SD/rno,fluxd,SDQ/rno
      endif
c     WRITE(2,*) XX(IN), fluxd*rno, fluxqd*rno, qflux
c     WRITE(3,*) XX(IN), fluxd, fluxh, fluxhi 
c     c      WRITE(3,*) XX(IN), fluxh,fluxhi,SCHQ/rno

      if (i_print .eq.1) then
         WRITE(3,*) SCHI/rno,SCHE/rno,SCHQ/rno,SD/rno,SDQ/rno
      endif

c     chi_conv = 2a/R
      chi_conv = 2/x_in(8)

c     sdi = (ne/ni)(Lni/Lne)De - z(nz/ni)(Lni/Lnz)Dz
      sdi = (1d0/x_in(27))*x_in(10)/x_in(11)*SD - 
     &      x_in(19)*x_in(20)/x_in(27)*x_in(22)/x_in(11)*SDQ

      x_out(1) = chi_conv*SD/rno
      x_out(2) = chi_conv*SCHE/rno
      x_out(3) = chi_conv*sdi/rno
      x_out(4) = chi_conv*SCHI/rno
      x_out(5) = chi_conv*SDQ/rno
      x_out(6) = chi_conv*SCHQ/rno
      x_out(7) = 0.0
      x_out(8) = 0.0

      x_out(1) = x_out(1)*x_in(10)
      x_out(2) = x_out(2)*x_in(12)
      x_out(3) = x_out(3)*x_in(11)*x_in(27)
      x_out(4) = x_out(4)*x_in(13)*x_in(27)*x_in(14)
      x_out(5) = x_out(5)*x_in(22)*x_in(20)
      x_out(6) = x_out(6)*x_in(13)*x_in(20)*x_in(14)

cc      WRITE(3,*) XX(IN), elflux,ionflux,qflux,sumflux
C  319 WRITE(*,330) SCHI*N*TE*2.*EI/(EN*R),SCHE*N*TE*2.*EE/(EN*R),
C     &SD*N*2./(EN*R),SCHQ,SDQ
      if (i_print .eq.1) WRITE(*,*)N,TE,R,EN,EE,EI
  330 FORMAT(/,2X,'CHIEFF=',G12.5,' CHEEFF=',G12.5,' DEFF=',G12.5,
     &' fluxd=', G12.5,' DQeff=',G12.5)
c     &' CHQEFF=',G12.5,' DQEFF=',G12.5)
C
c-----------------------------------------------------
c      WRITE(1,02305) SCHI,SCHE,XX(IN),IN
c      WRITE(2,02305) RFL,FLS,XX(IN),IN
02305 FORMAT(/2X,3G10.3,I5)
c-----------------------------------------------------
      RHOS=CS/OMEGAC
      DB=RHOS*CS
      RLT=R/(PR*LTH)
      CHINORM=1.5*SCHI*LNH*PR/(DB*RHOS)
c
      if (i_print .eq.1) WRITE(*,335) SHPE,CHIC,SCHEF,DEF,DELM
  335 FORMAT(' Electrom. transp. SHPE=',G11.3,' CHIC=',G11.3,
     &' SCHEF=',G11.3,' DEF=',G11.3,' DELM=',G11.3)
      if (i_print .eq.1) WRITE(*,336) CHINORM,RLT
  336 FORMAT(' CHINORM=',G11.3,' R/Lt=',G11.3)
      if (i_print .eq.1) WRITE(*,337) DB,RHOS
  337 FORMAT(' BOHM COEFF =',G11.3,' RHOS=',G11.3)
c
00150 CONTINUE

      end subroutine qfm_sub
  
