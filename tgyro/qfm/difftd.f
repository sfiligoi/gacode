C    DIFFTD calculates transport coefficients
C***************************************************************************
C     File DIFFTD.F
C     Assumes maximum 10 unstable roots
C
      SUBROUTINE DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DHQ)
C
      IMPLICIT NONE
      SAVE

      integer i_print

      INTEGER I,IR,IX,J,IW,LPRINTIN,NDIM,NEQ
      REAL*8 WR,WI,EN,ENI,ENH,EI,EE,TAUI,FT,TE,N,WIN,WDE
      REAL*8 XI(5),XE(5),XD(5),CHI(5),CHE(5),D(5)
      REAL*8 CHQ(5),DHQ(5),XQ(5),XDQ(5)
      REAL*8 U(5,100)
      REAL*8 DNI,DNE,WSQ,NN,DNIN,WRED,fre
      REAL*8 A,B,C,E,F,A1,DH,WR1,WI1,WR2,WI2
      REAL*8 SCHI,SCHE,SD,HP,SHP,NR,NI
      REAL*8 THRD,TVR,FTR,STR,FTRT,RFL,D1,XHH,DT,THRD2
      REAL*8 XIH,XEH,XDH,STF,PHS,KPC
      REAL*8 XDE,XDI,VEF,BTA,GAR,GAI,GBR,GBI,YDA
      REAL*8 HR,DIVGA,DIVGB,DIVA,DIVB,DEVGA,DEVGB,DEVA
      REAL*8 DEVB,SVB,WII,VEFN,BT1,EIH,EEH,EQH,IMP,TZ
      REAL*8 SCHQ,SDQ,KXQ,CHQEFF,DQEFF
      REAL*8 DTIQ,DTQQ,XQH,NQ,TQ,TIQ,NG,LNEH,LNHE
      REAL*8 BE,G,GI,SI,Z,TAUZ,EQ,ENQ,DTIMP,NIMP,H,ZFS
      REAL*8 KIQ,DNQ,DI,DE,DQ,NQR,NQI,K,T,DQT,DQN,TS
      REAL*8 WEXB,ROT,WROR,WIOR,WJR,WJI
      REAL*8 RVA,RVAID,RVAC,RVAC1,RVAC2,RVAC3
      REAL*8 RVAID1,RVAID2,RVAID3,RVB
      REAL*8 SVAID,SVAC,SVAC1,SVAC2,SVAC3
      REAL*8 SVAID1,SVAID2,SVAID3
      REAL*8 NIER,NIEI,RN
      REAL*8 DIVGA1,DIVGA2,DIVGA3,DIVA1,DIVA2,DIVA3
      REAL*8 DEVGA1,DEVGA2,DEVGA3,DEVA1,DEVA2,DEVA3
      REAL*8 DIVGAT,DIVAT,DEVGAT,DEVAT,SVAT,AT,CT,ET,HT,TT,DQNT
      REAL*8 CETAIN(32),ETE,ETI,ETQ,TAUH,AZ,FX,RFLT
      REAL*8 H1,H2,A2,A3,C1,C2,E1,E2,T1,T2,DQN1,DQN2
      REAL*8 HPT(5)
      REAL*8 PI,PE,PN,PQ,PNQ,PIQ,EA
      REAL*8 GCI,GCE,GD,GCQ,GNQ,GK
      REAL*8 LTI,LTE,LN,LTQ,LNQ
      REAL*8 ZVR(10,10),ZVI(10,10)
      COMPLEX*16 FIH,NEF,NRAT,ELN,TEF,FRAT,AV
      COMPLEX*16 NIE
      REAL*8 ELNR,ELNI,AINF,HPE,SHPE,CHIC,IMF,CEFT,SCHEF,GEI,DEFT,DEF
      REAL*8 betae,q,S,Cs,EM,EM1
      COMMON/LT/ LTI,LTE,LN,LTQ,LNQ
      COMPLEX*16 ZZ(10),RP(10),W,WJ,NE,BT2
      COMPLEX*16 GA,GB,GM,IU,HC
      COMMON/GRAD/ EN,ENH,EI,EE
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XHH,IX,IW
      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BE,G,SI,Z,TAUZ,ENQ,EQ,KIQ,KXQ,ZFS
      COMMON/WROT/ WEXB,ROT
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUH,AZ,FX,RFLT,
     1LPRINTIN,NDIM
      COMMON/TP/ EA,HPT
      COMMON/GRKVOT/ GK
      COMMON/ZV/ ZVR,ZVI
      COMMON/ZZ/ ZZ
      COMMON/NEQ/ NEQ
      COMMON/EM/ SHPE,CHIC,SCHEF,DEF
      COMMON/BETAE/ betae,q,S,Cs,EM
      COMMON/WDE/ WDE

      i_print = 0

C
C--------------------------------------------------------------
      EM1=EM
C      EM1=0.    !!! DEF not used if this line is active !
C-------------------------------------------------------------
C
      TVR=2.D0*THRD
      THRD2=THRD*THRD
      FTRT=FTR*TAUI
      IU=(0.D0,1.D0)
      SHP=0.D0
      CHIC=0.
      SHPE=0.
      SCHEF=0.
      DEF=0.
      ENI=1./EN
      EIH=EI-STR+FTR*ENH
      EEH=EE-STR+FTR*EN
      EQH=EQ-STR+FTR*ENQ
      LNEH=EN/ENH
      LNHE=ENH/EN
      G=1.-Z*BE-ZFS
      NG=MAX(G,1.D-10)
      GI=1./NG
      TZ=Z*TAUZ
      IF(TZ.LT.0.0001) GO TO 00005
      IMP=1./TZ
      KPC=1.
00005 CONTINUE
C
c      WRITE(*,00998) BE,G,SI,Z,TAUZ,ENQ,EQ,KIQ,KXQ,ZFS
c00998 FORMAT(' BQ=',G11.3,' G=',G11.3,' SI=',G11.3,' Z=',G11.3,/,
c     &' TAUZ=',G11.3,' ENQ=',G11.3,' EQ=',G11.3,' KIQ=',G11.3,/,
c     &' KXQ=',G11.3,' ZFS=',G11.3)
      DO 00010 I=1,5
      CHI(I)=0.D0
      CHE(I)=0.D0
      D(I)=0.D0
      CHQ(I)=0.D0
00010 DHQ(I)=0.D0
C
      SCHI=0.D0
      SCHE=0.D0
      SD=0.D0
      SCHQ=0.D0
      SDQ=0.D0
      XHH=0.D0
C
      DIVGA1=0.D0
      DIVGA2=0.D0
      DIVGA3=0.D0
      DIVGB=0.D0
      DIVA1=0.D0
      DIVA2=0.D0
      DIVA3=0.D0
      DIVB=0.D0
      DEVGA1=0.D0
      DEVGA2=0.D0
      DEVGA3=0.D0
      DEVGB=0.D0
      DEVA1=0.D0
      DEVA2=0.D0
      DEVA3=0.D0
      DEVB=0.D0
      RVA=0.D0
      RVAID=0.D0
      RVAID1=0.D0
      RVAID2=0.D0
      RVAID3=0.D0
      RVAC=0.D0
      RVAC1=0.D0
      RVAC2=0.D0
      RVAC3=0.D0
      SVAID=0.D0
      SVAID1=0.D0
      SVAID2=0.D0
      SVAID3=0.D0
      SVAC=0.D0
      SVAC1=0.D0
      SVAC2=0.D0
      SVAC3=0.D0
      SVB=0.D0
      PHS=0.D0
      RN=0.D0
C
      IF(NEQ.LT.7) VEF=0.D0
c
      DO 34 J=1,5
   34 HPT(J)=0.
C
      N=U(3,IX)
      TE=U(2,IX)
      TQ=U(4,IX)
      NQ=U(5,IX)
C
C  IR IS THE NUMBER OF UNSTABLE MODES
C  A SUMMATION IS PERFORMED OVER THE DIFFUSION
C  DUE TO ALL UNSTABLE MODES
C
c      WRITE(*,20011) EA
20011 FORMAT(/,' DIFFTD  EA=',G11.3)
      fre=0.25
c      fre=0.5
      IF(IR.EQ.0) GOTO 00110
c  ---------------------------------------------------
c  Main Loop *******************
c
      DO 01100 J=1,IR
      WJ=RP(J)
      WJR=DREAL(WJ)
      WJI=DIMAG(WJ)
      IF(WJI.LT.0.001) GO TO 01100
c
      IF(IW.NE.1) GO TO 10022
      WROR=WJR*WDE
      WIOR=WJI*WDE
      if (i_print .eq. 1) WRITE(*,10021) WROR,WIOR,J
10021 FORMAT(' Orig W, WR=',G11.3,'/s WI=', G11.3,'/s  J=',I5)
10022 CONTINUE
c
      W=RP(J)-IU*ROT*DABS(WEXB)
      WR=DREAL(W)
      IF(WR.LE.0.) GO TO 10011
      WRED=(DIMAG(RP(J)))**2-fre*ROT*WEXB**2
      IF(WRED.LT.0.D0) WRED=0.D0
      W=WR+IU*DSQRT(WRED)
C      W=RP(J)    !!!!!!!  NOTE TEMPORARY CHANGE  NO STAB OF TE MODE
10011 CONTINUE
      WI=DIMAG(W)
      IF(WI.GT.0.001) GO TO 11 
      WI=0.001
      GO TO 01100
   11 CONTINUE
      W=DCMPLX(WR,WI)
      GM=1.D0+EE*ENI/(W-1.D0+IU*VEF)
      BT1=BTA-2.5D0
      BT2=BTA-2.5D0*GM
      HC=W-FTR+TVR*BT1
      NE=W*W-2.D0*FTR*W+FTR+IU*VEF*HC
      NIE=W*W-2.D0*FTR*W+FTR
      NIER=DREAL(NIE)
      NIEI=DIMAG(NIE)
      WSQ=WR*WR+WI*WI
      WII=1.D0/WI
C
C
C   ******   COLLISION PARAMETERS  *******
C
      GA=W-FTR+TVR*BT2
      GB=(W-FTR)/(W-1.D0+IU*VEF)
      GAR=DREAL(GA)
      GAI=DIMAG(GA)
      GBR=DREAL(GB)
      GBI=DIMAG(GB)
      HR=WR-FTR+TVR*BT1
      XDE=WSQ-FTR*WR
      XDI=WSQ+FTR*TAUI*WR
      YDA=WR*(1.D0-EN)+EE-STR+FTR*EN
c Linear trapped electron density response dn?n = FX ephi/Te where **
c **  FX = ENI*(YDA+IU*WI*(1-EN)+IU*VEF*(EN*GA+EE*GB))/NE **
C ** where NE=NER+IU*NEI; NER=NIER-WI*VEF, NEI=NIEI+VEF*HR ****
C   ***************************************
C
      IF(IW.NE.1) GOTO 00021
c      WR1=DABS(EN)*WR
c      WI1=DABS(EN)*WI
      WR2=WDE*WR
      WI2=WDE*WI
      WR1=WR
      WI1=WI
      if (i_print .eq. 1) WRITE(*,00020) WR1,WI1,J
00020 FORMAT(2X,'WR=',D11.5,' WI=',D11.5,' J=',I5)
      if (i_print .eq. 1) WRITE(*,10020) WR2,WI2,J
10020 FORMAT('  WR in sec**-1 ',G11.3,' WI in sec**-1 ',G11.3,' J=',I5)
00021 CONTINUE
C
      NR=DREAL(NE)
      NI=DIMAG(NE)
      NN=(NR)**2+(NI)**2
C
      IF(VEF.EQ.0.D0) GO TO 25
c
c    We write FX = KK*(KR + IU*KI) where KK=1/(NN*EN),
c    KR = RVA + EE*RVB and KI = SVA + EE*SVB
c    We divide into ideal and collisional parts as:
c    RVA = RVAID + RVAC, SVA = SVAID + SVAC 
c
      RVAID=NIER*YDA+NIEI*W*(1.D0-EN)
      SVAID=NIER*WI*(1.D0-EN)-NIEI*YDA
c
      RVAC=VEF*(EN*(NI*GAR-NR*GAI)-WI*YDA+HR*WI*(1.D0-EN))
      RVB=VEF*(NI*GBR-NR*GBI)
c
      SVAC=VEF*(EN*(NR*GAR+NI*GAI)-(1.D0-EN)*WI**2-HR*YDA)
      SVB=VEF*(NR*GBR+NI*GBI)
c
c     These parts are now divided into diagonal, off diagonal and convective 
c     parts as e.g.  RVA = RVA1 + EN*RVA2 + EE*RVA3 etc ...
c
c ** The following parts, due to the ideal part of the density responce
c    enter only for electron thermal transport -----------------------
c
      RVAID1=NIER*(WR-STR)+WI*NIEI
      RVAID2=NIER*(FTR-WR)-WI*NIEI
      RVAID3=NIER
c
      SVAID1=WI*NIER-NIEI*(WR-STR)
      SVAID2=-WI*NIER+NIEI*(WR-FTR)
      SVAID3=-NIEI
c----------------------------------------------------------------------
c
      RVAC1=VEF*WI*(HR-WR+STR)
      RVAC2=VEF*(GAR*NI-GAI*NR+WI*(WR-FTR-HR))
      RVAC3=-VEF*WI
c
      SVAC1=VEF*(HR*(STR-WR)-WI*WI)
      SVAC2=VEF*(GAR*NR+GAI*NI+WI*WI+HR*(WR-FTR))
      SVAC3=-VEF*HR
c
c     The ideal parts RVAID etc .. will only enter into the electron thermal 
c     transport. For the particle transport we need only SVAC and SVB
c
c  Ion thermal conductivity
c
      DIVGA1=FTRT*RVAC1
      DIVGA2=FTRT*RVAC2
      DIVGA3=FTRT*RVAC3
c
      DIVA1=XDI*SVAC1
      DIVA2=XDI*SVAC2
      DIVA3=XDI*SVAC3
c
      DIVGB=FTRT*RVB
      DIVB=XDI*SVB
c
c   Electron thermal conductivity
c
      DEVGA1=-FTR*RVAC1-BT1*VEF*(SVAID1+SVAC1)
      DEVGA2=-FTR*RVAC2-BT1*VEF*(SVAID2+SVAC2)
      DEVGA3=-FTR*RVAC3-BT1*VEF*(SVAID3+SVAC3)
c
      DEVA1=XDE*SVAC1-BT1*VEF*(WR-FTR)*(RVAID1+RVAC1)
      DEVA2=XDE*SVAC2-BT1*VEF*(WR-FTR)*(RVAID2+RVAC2)
      DEVA3=XDE*SVAC3-BT1*VEF*(WR-FTR)*(RVAID3+RVAC3)
c
      DEVGB=-(FTR*RVB+BT1*VEF*SVB)
      DEVB=XDE*SVB-BT1*VEF*(WR-FTR)*RVB
c
      PHS=(WR-FTR)*DREAL(BT2)+WI*DIMAG(BT2) !! independent of density resp.
c
      RN=1.D0/NN

      DIVGA1=DIVGA1*RN
      DIVGA2=DIVGA2*RN
      DIVGA3=DIVGA3*RN
      DIVGB=DIVGB*RN
      DEVGA1=DEVGA1*RN
      DEVGA2=DEVGA2*RN
      DEVGA3=DEVGA3*RN
      DEVGB=DEVGB*RN
      DIVA1=DIVA1*RN
      DIVA2=DIVA2*RN
      DIVA3=DIVA3*RN
      DIVB=DIVB*RN
      DEVA1=DEVA1*RN
      DEVA2=DEVA2*RN
      DEVA3=DEVA3*RN
      DEVB=DEVB*RN
      SVAC1=SVAC1*RN
      SVAC2=SVAC2*RN
      SVAC3=SVAC3*RN
      SVB=SVB*RN
   25 CONTINUE
C
      STF=STR-FTR*EN
      A1=WSQ*(EN-1.D0)+2.D0*WR*STF
      A=WSQ*(A1+FTR*(STR*EN-11.D0*THRD-TAUI*(1.D0-
     1FTR*EN)))+FTR*FTR*TAUI*(2.D0*WR*(1.D0-EN)-STF)
      A=A/NN
      B=(WSQ*(2.D0*(WR-FTR)+FTR*TAUI)-
     1FTR*FTR*TAUI)/NN
      C=WSQ*(A1+TVR*FTR*(EN-4.D0))+
     1FTR*FTR*(2.D0*WR*(EN-1.D0)+STF)
      C=C/NN
      DH=(WSQ*(2.D0*WR-5.D0)+FTR*FTR)/NN
      E=WSQ*(1.D0-EN)-2.D0*WR*STF+FTR*(11.D0*THRD
     1-STR*EN)
      E=E/NN
      F=2.D0*(-WR+FTR)/NN
C
      NQR=WR**2-WI**2+2.D0*FTR*IMP*WR+FTR*IMP*IMP
      NQI=2.D0*WI*(WR+FTR*IMP)
      NIMP=NQR**2+NQI**2
C
C   ****   SPLITTING IN  EN  AND EE  ************
C
      A2=(WSQ*WSQ+FTR*((STR+FTR*TAUI-2.D0*WR)*WSQ+FTR*TAUI*(FTR-
     &2.D0*WR)))/NN
      A3=(WSQ*(-WSQ+2.D0*STR*WR-FTR*(11.D0*THRD+TAUI))+FTR*FTR
     &*TAUI*(2.D0*WR-STR))/NN
      C1=(NN-WSQ*WSQ+14.D0*THRD*WSQ*WR-40.D0*THRD2*WSQ
     &-50.D0*THRD2*WR+175.D0/27.D0)/NN
      C2=(WSQ*WSQ-10.D0*THRD*WSQ*WR+10.D0*THRD2*WSQ
     &+50.D0*THRD2*WR-125.D0/27.D0)/NN
      E1=(WSQ+11.D0*THRD*FTR-2.D0*WR*STR)/NN
      E2=-(-2.D0*WR*FTR+(WSQ+STR*FTR))/NN
      H1=(WSQ*(-WSQ-2.D0*IMP*WR*STR
     &+FTR*IMP*IMP*(-11.D0*THRD)+FTR*TAUI*IMP)
     &+FTR*FTR*TAUI*IMP*IMP*(2.D0*WR+STR*IMP))/NIMP
      H2=(WSQ*(WSQ+2.D0*IMP*WR*FTR+FTR*IMP*IMP*STR-FTR*TAUI*IMP*FTR)
     &-FTR*FTR*TAUI*IMP*IMP*(2.D0*WR+FTR*IMP))/NIMP
      T1=KXQ*(WSQ*(-WSQ-2.D0*IMP*WR*STR-FTR*IMP*IMP*8.D0*THRD)
     &+FTR*FTR*IMP**3*(2.D0*WR+IMP*STR))/NIMP
      T2=KXQ*(WSQ*(WSQ+2.D0*IMP*WR*FTR+FTR*IMP*IMP*TVR)
     &-FTR*FTR*IMP**3*(2.D0*WR+IMP*FTR))/NIMP
      DQN1=(-NQR+2.D0*(WR+IMP*STR)*(WR+FTR*IMP))/NIMP
      DQN2=(NQR-2.D0*(WR+IMP*FTR)*(WR+FTR*IMP))/NIMP
C
      DIVGAT=DIVGA1+EN*DIVGA2+EE*DIVGA3
      DIVAT=DIVA1+EN*DIVA2+EE*DIVA3
      DEVGAT=DEVGA1+EN*DEVGA2+EE*DEVGA3
      DEVAT=DEVA1+EN*DEVA2+EE*DEVA3
      SVAT=SVAC1+EN*SVAC2+EE*SVAC3
      AT=EN*A2+A3
      CT=C1-1.D0+EN*C2
      ET=E1+EN*E2
      HT=H1+ENQ*H2
      TT=T1+ENQ*T2
      DQNT=DQN1+EN*DQN2
C
C **** IMPURITIES *****
C
      NIMP=(WR*(WR+2.*FTR*IMP)-WI*WI+FTR*IMP*IMP)**2
     1+4.D0*WI*WI*(WR+FTR*IMP)**2
      DTIMP=(WSQ*(WSQ*(ENQ-1.D0)+2.*IMP*WR*EQH+FTR*IMP
     1*IMP*(2.D0*EQ-11.D0/3.D0+STR*ENQ)+FTR*TAUI*IMP*(1.D0+
     1EQ-FTR*ENQ))+FTR*TAUI*(2.D0*FTR*IMP*IMP*WR*(1.D0-ENQ)
     1-FTR*IMP**3*EQH))/NIMP
C  *************
C
      DNI=(WR+FTR*TAUI)**2+WI*WI
      DNE=(WR-FTR)**2+WI*WI
      DNQ=(WR+FTR*IMP)**2+WI*WI
      DI=EI*RFL*DNI
      DE=EE*RFL*DNE
      DQ=EQ*RFL*DNQ
C
C **** IMPURITIES *****
C
      NIMP=NQR**2+NQI**2
C
      H=(WSQ*(WSQ*(ENQ-1.D0)-2.D0*IMP*WR*(STR-FTR*ENQ)
     1+FTR*IMP*IMP*(-11.D0*THRD+STR*ENQ)+FTR*TAUI*IMP*(1.D0-FTR
     1*ENQ))+FTR*FTR*TAUI*IMP*IMP*(2.D0*WR*(1.D0-ENQ)+(STR-FTR
     1*ENQ)*IMP))/NIMP
      K=IMP*(FTR*FTR*TAUI*IMP*IMP-WSQ*(2.D0*WR
     1+FTR*(2.D0*IMP+TAUI)))/NIMP
C
C  *************
      T=(WSQ*(WSQ*(ENQ-1.D0)-2.D0*IMP*WR*(STR-FTR*ENQ)
     1+FTR*IMP*IMP*(-8.D0*THRD+TVR*ENQ))+FTR*FTR*(IMP)**3
     1*(2.D0*WR*(1.D0-ENQ)+IMP*(STR-FTR*ENQ)))/NIMP
      TS=IMP*(FTR*FTR*(IMP)**3-WSQ*(2.D0*WR+5.D0*IMP))/NIMP
C
      T=KXQ*T
      TS=KXQ*TS
      DQN=(ENQ-1.D0)*NQR+2.D0*((1.D0-ENQ)*WR+IMP*(STR-FTR*ENQ)
     1)*(WR+FTR*IMP)
      DQT=2.D0*IMP*(WR+FTR*IMP)
      DQN=DQN/NIMP
      DQT=DQT/NIMP
C
      DTIQ=H-EQ*K
      DTQQ=T-EQ*TS
C
      XDH=D1*(TE**1.5D0)*WI**3/RFL
      XIH=XDH/DNI
      XEH=XDH/DNE
      XQH=XDH/DNQ
C
C
      XI(1)=XIH
      XI(2)=TVR*FT*XIH*TAUI*GI*(B-DIVGA3-DIVGB-WII*(DIVA3+DIVB))
      XI(3)=-TVR*XIH*(TE/N)*TAUI*(LNEH+FT*GI*(A3+DIVGA1+DIVA1*WII))
      XI(4)=-XIH*TVR*BE*Z*GI*K*TAUZ*TAUI
      XI(5)=TVR*XIH*(TE/N)*TAUI*Z*H1*GI
C
      PI=FT*TVR*XIH*GI*(A2+DIVGA2+WII*DIVA2)*EA
      PIQ=-TVR*XIH*Z*H2*BE*GI*EA
C
C
      XE(1)=0.D0
      XE(2)=FT*XEH*(1.D0+TVR*(DH-DEVGB-DEVGA3-WII*(DEVB+DEVA3)))
      XE(3)=-TVR*FT*XEH*(C1+DEVGA1+WII*DEVA1)*TE/N
      XE(4)=0.D0
      XE(5)=0.D0
C
      PE=FT*TVR*XEH*(C2+DEVGA2+WII*(DEVA2+VEF*PHS))*EA
C
C
      XD(1)=0.D0
      XD(2)=-FT*XDH*(N/TE)*(F+WII*(SVB+SVAC3))
      XD(3)=FT*XDH*(E1-WII*SVAC1)
      XD(4)=0.D0
      XD(5)=0.D0
C
      PN=-FT*XDH*(E2-WII*SVAC2)*EA
C
C
      XQ(1)=0.D0
      XQ(2)=0.D0
      XQ(3)=0.D0
      XQ(4)=(1.D0+TVR*TS)*XQH
      XQ(5)=-TVR*TQ*(1.D0+T1)*XQH/NQ
C
      PQ=TVR*XQH*T2*EA
C
C
      XDQ(1)=0.D0
      XDQ(2)=0.D0
      XDQ(3)=0.D0
      XDQ(4)=-XDH*DQT*NQ/TQ
      XDQ(5)=XDH*DQN1
C
      PNQ=-XDH*DQN2*EA
C
      TIQ=1.D0/(TAUZ*KIQ)
C
      HP=XIH*GI*TAUI*TVR*FTR*(1.D0-FT)*EA
      GCI=XI(1)+(LNHE*(XI(2)*EE+XI(3)*N/TE)+XI(4)*EQ*TIQ
     1+XI(5)*NQ/(TE*KIQ))/(EI*TAUI)-(HP+PI+PIQ)*GK*2.D0*LTI
      GCE=XE(2)+XE(3)*N/(TE*EE)-PE*GK*2.D0*LTE
      GD=XD(2)*TE*EE/N+XD(3)-PN*GK*2.D0*LN
      GCQ=XQ(4)+XQ(5)*NQ/(TQ*EQ)-PQ*GK*2.D0*LTQ
      GNQ=XDQ(4)*TQ*EQ/NQ+XDQ(5)-PNQ*GK*2.D0*LNQ
C
C
      SHP=SHP+HP
C
      HPT(1)=HPT(1)+HP+PI+PIQ
      HPT(2)=HPT(2)+PE
      HPT(3)=HPT(3)+PN
      HPT(4)=HPT(4)+PQ
      HPT(5)=HPT(5)+PNQ
C
      SCHI=SCHI+GCI
      SCHE=SCHE+GCE
      SD=SD+GD
      SCHQ=SCHQ+GCQ
      SDQ=SDQ+GNQ
C
C
      CHQEFF=D1*TE**1.5*WI**3*(EQ-TVR-TVR*DTQQ)/DQ
      DQEFF=XDH*(DQN-EQ*DQT)
00222 FORMAT(2X,'CHQEFF=',G11.4,' DQEFF=',G11.4,' XIH=',G11.3,
     &' XEH=',G11.3)
C
      IF(IW.NE.1) GO TO 00093
       if (i_print .eq. 1) WRITE(*,00222) CHQEFF,DQEFF,XIH,XEH
       if (i_print .eq. 1) WRITE(*,00099) SCHI,SCHE,SD,SCHQ,SDQ,J
00099 FORMAT(2X,'SCHI=',F11.5,' SCHE=',F11.5,' SD=',
     1F11.5,' SCHQ=',F11.5,' SDQ=',F11.5,' J=',I5)
       if (i_print .eq. 1) WRITE(*,00223) XDH,DNE,WI
00223 FORMAT(' XDH=',G11.3,' DNE=',G11.3,' WI=',G11.3)
C      WRITE(*,00229) XEH,DH,DEVGB,DEVGA3,WII,DEVB,DEVA3,XE(2)
00229 FORMAT(' XEH=',G11.3,' DH=',G11.3,' DEVGB=',G11.3,' DEVGA3=',G11.3,
     &/,' WII=',G11.3,' DEVB=',G11.3,' DEVA3=',G11.3,' XE(2)=',G11.3)
00093 CONTINUE
C
      XHH=XHH+XIH
C
      DO 00100 I=1,5
      CHI(I)=CHI(I)+XI(I)
      CHE(I)=CHE(I)+XE(I)
      D(I)=D(I)+XD(I)
      CHQ(I)=CHQ(I)+XQ(I)
      DHQ(I)=DHQ(I)+XDQ(I)
00100 CONTINUE
C
C
01100 CONTINUE
c-----------------------------------------------  MAIN LOOP
C
C--------------------------------------------------------------
      IF(NEQ.LE.8) GO TO 00110
C 
      IF(EM.EQ.0.AND.EM1.EQ.0.) GO TO 00110
C  IF electromagnetic effects or collisions on free electrons are included
C  the transport coefficients are corrected for this.
C------------------------------------------------------------
      DT=D1*TE**1.5D0
      SHPE=0.D0
      SCHEF=0.D0
      DEF=0.D0
      DO 00200 J=1,NEQ
      WR=DREAL(ZZ(J))
      IF(DIMAG(ZZ(J)).LE.0.01D0) GO TO 00200
      WI=DIMAG(ZZ(J))-ROT*DABS(WEXB)
      W=WR+IU*WI
      IF(WR.LE.0.D0) GO TO 10012
      WRED=(DIMAG(ZZ(J)))**2-fre*ROT*WEXB**2
      IF(WRED.LT.0.D0) WRED=0.D0
      W=WR+IU*DSQRT(WRED)
C      W=ZZ(J)  !!!!  NOTE TEMPORARY CHANGE NO STAB OF TE MODE
10012 CONTINUE
      WI=DIMAG(W)
      WIN=DIMAG(ZZ(J))
      IF(WIN.LE.1.D-3) WIN=1.D-3
      IF(WI.LT.1.D-3) GO TO 00200
c --- contr. to chii from em free electr. ----
c
      DNI=(WR+FTR*TAUI)**2+WI*WI
      DNIN=(WR+FTR*TAUI)**2+WIN*WIN
      XDH=DT*WI**3/RFL
      XIH=XDH/DNI
      FIH=DCMPLX(ZVR(1,J),ZVI(1,J))
      IF(NEQ.EQ.11) GO TO 197
C----------------------------------------------
      IF(NEQ.EQ.9) GO TO 00095
c-- Here NEF for disp10 is defined ---
      NEF=DCMPLX(ZVR(4,J),ZVI(4,J))
      GO TO 00097
C----------------------------------------
00095 CONTINUE
C-- Here NEF for disp 9 is defined --
      AV=DCMPLX(ZVR(8,J),ZVI(8,J))
      NEF=FIH-(ZZ(J)-ENI)*AV/KPC
C--------------------------------------------
      GO TO 00097
  197 CONTINUE
      AV=DCMPLX(ZVR(9,J),ZVI(9,J))
      NEF=FIH-(ZZ(J)-ENI)*AV/KPC
      TEF=EE*ENI*AV/KPC
00097 CONTINUE
      IF(CDABS(FIH).LT.0.0001) FIH=(0.0001,0.)
      NRAT=NEF/FIH
      ELN=NRAT-1.D0
      ELNR=DREAL(ELN)
      ELNI=DIMAG(ELN)
      AINF=TVR*(FTR*TAUI*ELNR+ELNI*(WI*WI+WR*(WR+FTR*TAUI))/WIN)
      HPE=XIH*GI*(1.D0-FT)*EA*AINF
      IF(IW.NE.1) GO TO 20035
C      WRITE(*,20021) ELNR,ELNI,WR,WI,WIN,AINF,HPE
20021 FORMAT(' ELNR=',G11.3,' ELNI=',G11.3,' WR=',G11.3,' WI=',G11.3,
     &' WIN=',G11.3,' AINF=',G11.3,' HPE=',G11.3)
20035 CONTINUE
      SHPE=SHPE+HPE
c
c ****  Free electron heat flux *********************
c
C-----------------------------------------------------
      IF(NEQ.EQ.11) GO TO 10099
      IF(NEQ.EQ.9) GO TO 10098
c--- Here TEF for disp10 is defined ---
      TEF=DCMPLX(ZVR(6,J),ZVI(6,J))
      GO TO 10099
10098 CONTINUE
c-- Here TEF for disp9 is defined ---
      TEF=EE*ENI*AV/KPC
C-------------------------------------------------------------
10099 CONTINUE
      FRAT=TEF/FIH
      IMF=-DIMAG(FRAT)*DNIN/DNI
      CEFT=(1.-FT)*IMF*DT*ETE*WI**3/(RFL*WIN)
      SCHEF=SCHEF+CEFT
      IF(IW.NE.1) GO TO 20036
C      WRITE(*,20022) FRAT,IMF,CEFT,SCHEF
20022 FORMAT(' FRAT=',2G11.3,' IMF=',G11.3,' CEFT=',G11.3,
     &' SCHEF=',G11.3)
20036 CONTINUE
C**********************************************************
c ----Free electron particle flux -----------
      GEI=-DIMAG(NRAT)/WI
      DEFT=(1.D0-FT)*GEI*EN*XDH
      DEF=DEF+DEFT
c -------
00200 CONTINUE
c
      CHIC=-SHPE*GK*2.D0*LTI
      HPT(1)=HPT(1)+EM*SHPE
      CHE(2)=CHE(2)+EM*SCHEF
      D(3)=D(3)+EM1*DEF
c
C-------------------------------------------------------
00110 CONTINUE
c
c      WRITE(*,21001) HPT(1),HPT(2),HPT(3),HPT(4),SHPE
21001 FORMAT(' HP1',G11.3,' HP2',G11.3,' HP3',G11.3,' HP4',G11.3,
     &' SHPE=',G11.3)
C
      RETURN
      END
