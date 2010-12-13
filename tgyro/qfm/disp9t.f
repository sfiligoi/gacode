C***********************************************************************
C  Linear solver of dispersion equation in matrix form
C**********************************************************************
      subroutine disp9t(neq,ZZ)
c
c  ***************************************************
c  THIS ROUTINE IS A MODIFICATION OF THE LINEAR PART OF ETAWN6
c  WRITTEN BY GLENN BATEMAN. THE MODIFICATIONS CONSIST OF THE INCLUSION
c  OF IMPURITIES IN THE DILUTION APPROX. IN THE SYSTEM WITH
c  4 EQUATIONS AND THE INCLUSION OF COLLISIONS ON TRAPPED
c  ELECTRONS IN THE FULL SYSTEM. THIS SYSTEM THEN CONSISTS
c  OF 7 EQUATIONS. WHEN PARALLEL ION MOTION IS INCLUDED THERE
c  ARE 8 EQUATIONS WITH COLLISIONS. IN disp10 ALSO ELECTROMAGNETIC 
c  EFFECTS ARE INCLUDED.  
c  MOST PARAMETERS ARE TRANSFERRED THROUGH COMMON BLOCKS LIKE GRAD,
c  IMP AND ETAWN6. NOTE THE INVERSE DEFINITION OF ZTAUH AND ZTAUZ !
c  ********************************************************
      implicit none
c      SAVE
c
      INTEGER idp
c------------------------------------------------------------------
c  Note  that idp gives the dimensions of zvr and zvi which also have to be
c  defined in the main programme and in difftd.f. These declarations thus
c  have to be  the same.
      PARAMETER ( idp = 10)
c-----------------------------------------------------------------------
c
      LOGICAL inital, lmatv
      data inital /.true./, lmatv /.true./
c
      REAL*8  cetain(32)
     &  , omega(idp), gamma(idp)
      COMPLEX*16 ZZ(*), ZZS(10)
c
      REAL*8 epsnhin, epsnzin, epstein, epsthin, epstzin, tauhin, tauzin
     & , fnzin, czin, azin, ftrapein, epsnin
     & , ekyrhoin, g,   etai, etae, rb2, rb3, W0I, W0Z, sigp, sigt
     & , zb, si, eq, kiq, kxq, bt, bt1, vef, tvr, ftr, ftrt, rbou
c
      REAL*8 KAPPA, RAV, GAV,ALA,SH,SH2,WZR,WZI,WZIMAX
      REAL*8 H1,XH,EN,EI,EE,TAU,TAUI,FL,FT,GM,BTA,XT,R,HQR,HQI,WM
      COMPLEX*16 ALPC,ALPK,ALPHA,HQ,WZ,WZ1,WZ2,IU,H2,E,E1,WZJ(100),WZH
      COMPLEX*16 WZP,WZPP,DWN,DEW,DEWP,D2EW,W1,W2,AM1,BM1,DR
      COMPLEX*16 DWN1,WZA,WS,WSF,WSFF,WZO,WZE
      REAL*8 DW1,DW2,CTEST
c
      SAVE WZJ,idim,tvr,ftr,em1,A,zepsilon,zepsmach
      INTEGER lprintin, lprint, neq, idim, ndim 
     & , ieq, j1, j2, j, iret, IK, IST, IM, ISB,IMX
c
c ndim  = first dimension of the 2-D array difthi
c           and the maximum number of unstable modes allowed
c ieq   = number of equations
c
      REAL*8 zamr(idp,idp),zami(idp,idp),zbmr(idp,idp),zbmi(idp,idp)
     &  ,zamrt(idp,idp), zbmrt(idp,idp)
     &  ,zalfr(idp),zalfi(idp),zbeta(idp),zvr(idp,idp),zvi(idp,idp),ztol
c
      INTEGER iter(idp), ifail
c
c zamr(i,j) = matrix A
c zbmr(i,j) = matrix B
c   Note that the eigenvalues are
c omega(j) = ( zalfr(j) + i zalfi(j) ) / zbeta(j)
c where beta(j) will be 0.0 in the case of an infinite eigenvalue
c zvr(j) = eigenvector real part
c zvi(j) = eigenvector imag part
c
      REAL*8 wr,wi,H,fft,fzft
c
      REAL*8 zepsilon, zepsmach, zepsqrt
     & , zetah, zetaz, zetae 
     & , zepsnh, zepste, zepsth
     & , ztauh, ztauz, zep2nh, zep2nz, zep2ne, zft
     & , zimp, zfnz, zmass, zflh, zflz,zfs,A,em,em1
      REAL*8 q,S,Cs,alp,k1,k2,kps,betae,eni,alf,kpc
      INTEGER ITC,ITL,ITS,ITERA,IMET,ISEARCH,SEARCHMODE
      REAL*8 TOL
c
c
      COMMON/GRAD/ epsnin,epsnhin,etai,etae
      COMMON/IMP/ fnzin,g,si,czin,tauzin,epsnzin,eq,kiq,kxq,zfs
      COMMON/ETAWN6/ cetain,epstein,
     &epsthin,epstzin,tauhin,azin,
     &ftrapein,ekyrhoin,lprintin,ndim
      COMMON/IRET/ iret
      COMMON/COLL/ bt,vef
      COMMON/BETAE/ betae, q, S, Cs, em
      COMMON/TEST/ alp,alf,kpc,kps
      COMMON/HQ/ H1,HQ,GAV,WZ,WZP,WZIMAX
      COMMON/SHAFRS/ ALA
      COMMON/KAPPA/ KAPPA,RAV
      COMMON/IK/ IK,IST,ITC,ITL,ITS,ITERA,TOL
      COMMON/ZV/ zvr,zvi
C      COMMON/FZ/ fft,faft,ftt,ftf,fzf,fh
      COMMON/ISB/ ISB
      COMMON/SEARCHMODE/ SEARCHMODE
c
c ...  STRUCTURE OF MATRIX EQUATION ...
c
c ...  omega*zbmr(i,j) = zamr(i,j)
c
c    variables i=1,6: efi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
c    variables j=1,6 same as for i
c
c  ...................................................
c
c  ---------------------------------------------------
c  Iteration control variables used only when neq=9.
c
c   ITC=1  means  iterate ITC=2 means use only average
c   ITS=1  means  iterations converged
c   ITL    means  maximum number of iterations
c   IST=1 means  start iterations from analytical approx.
c   ITER   means  number of iterations performed
c   TOL    means  relative error allowed for convergence
c   IF IST.NE.1 the iterations will start from whatever value stored in 
c   WZJ(IK). This will during the simulation be the root from the previous time
c   step with the largest growthrate corresponding to propagation in the ion 
c   drift direction.
c   It is necessary that IST=1 in the first time step
c  --------------------------------------------------------
c
c..initialize variables
c
      if ( inital ) then
c
        idim = idp
c
      if(neq.LE.idp) go to 40001
      write(*,40002)
40002 format(' The dimension idp is not sufficient') 
      go to 00085
40001 continue
        tvr = 2./3.
        ftr = 5./3.
        ftrt = ftr/tauhin
        em1=em
c        em1=1.D0
c
        A = azin
        zepsilon = 1.e-4
c
        zepsmach = 0.5
  2     if ( 0.5 * zepsmach + 1.0 .gt. 1.0 ) then
          zepsmach = 0.5 * zepsmach
          go to 2
        endif
c
        zepsqrt = sqrt ( zepsmach )
c
        inital = .false.
c
      endif
c
      lprint = lprintin
c
      ftrt=ftr/tauhin
      A = azin
      em1 = em
      ieq = max ( 2, neq )

c
      IU=(0.,1.)
c
      ftrt = ftr/tauhin
      bt1 = bt - 2.5
c
      iret = 0
c..print header
c
      if ( lprint .gt. 2 ) then
c
        write (6,*)
        write (6,*)
     & 'Weiland-Nordman eigenvalue equations, subroutine etawn6'
        write (6,*) '(all frequencies normalized by omega_{De})'
        write (6,*) '(all diffusivities normalized by '
     &    ,'omega_{De} / k_y^2'
c
      endif
c
c..check validity of input data
c
c      if ( neq .lt. 2 ) call abortb (6
c     & ,' neq .lt. 2 in sbrtn etawn6')
c
c      if ( abs(epstein) .lt. zepsqrt ) call abortb (6
c     & ,' abs(epstein) .lt. zepsqrt in sbrtn etawn6')
c
c      if ( abs(epsthin) .lt. zepsqrt ) call abortb (6
c     & ,' abs(epsthin) .lt. zepsqrt in sbrtn etawn6')
c
c      if ( neq .gt. 4 ) then
c
c        if ( abs(epstzin) .lt. zepsqrt ) call abortb (6
c     &   ,' abs(epstzin) .lt. zepsqrt in sbrtn etawn6')
c
c        if ( abs(epsnzin) .lt. zepsqrt ) call abortb (6
c     &   ,' abs(epsnzin) .lt. zepsqrt in sbrtn etawn6')
c
c        if ( czin .lt. 1.0 ) call abortb (6
c     &   ,' czin .lt. 1.0 in sbrtn etawn6')
c
c      endif
c
      do 10 j1=1,ndim
        omega(j1) = 0.0
        gamma(j1) = 0.0
  10  continue
c
      zepsth = epsthin
      zepsnh = epsnhin 
      zepste = epstein 
c
c..compute the rest of the dimensionless variables needed
c
      zetah  = zepsnh / zepsth
*USCH      zetae  = zepsne / zepste
c
c  *******  NOTE THE INVERSE DEFINITION OF ZTAUH ! ******
      ztauh  =1./tauhin
c
      zep2nh=epsnhin
      zep2ne=epsnin
      zetae = zep2ne/zepste
      zft    = ftrapein
      zflh   = ekyrhoin**2
      eni = 1./zep2ne
c
      zimp   = czin
      zmass  = azin
c
        zfnz   = fnzin * zimp
c        zetaz  = zepsnz / zepstz
        zetaz=epsnzin/epstzin
c
c  ******  NOTE THE INVERSE DEFINITION OF ZTAUZ ! ******
        ztauz = 1./(czin*tauzin)
        zep2nz=epsnzin
c        zep2nz = 2.0 * zepsnz
        zflz   = zmass * zflh / zimp**2
        rbou=1.d0
c        rbou=0.25d0+S*2.d0/3.d0
c       !!!Bounce average of Omega_De !!!
c        rb2=2./3.+5.*S/9.
        rb2=1.0
c       !!!Average of omega_Di!!!
c        rb3=1.+0.79*(S**2)
        rb3=1.0
c       !!!Average of FLR!!!
        W0Z=0.d0
        W0I=0.d0
        sigp=0.D0
        sigt=0.D0
c..diagnostic output
c
      if ( lprint .gt. 6 ) then
        write (6,*)
        write (6,*) '--------------------------------------'
        write (6,*)
        write (6,*)
        write (6,*) zetah,' = zetah'
        write (6,*) zetaz,' = zetaz'
        write (6,*) zetae,' = zetae'
        write (6,*) ztauh,' = ztauh'
        write (6,*) ztauz,' = ztauz'
        write (6,*) zep2nh,' = zep2nh'
        write (6,*) zep2nz,' = zep2nz'
        write (6,*) zep2ne,' = zep2ne'
        write (6,*) zft,' = zft'
        write (6,*) zimp,' = zimp'
        write (6,*) zmass,' = zmass'
        write (6,*) zfnz,' = zfnz'
        write (6,*) zflh,' = zflh'
        write (6,*) zflz,' = zflz'
        write (6,*) zepsqrt,' = zepsqrt'
        write (6,*) zepsmach,' = zepsmach'
        endif
c------------------------------------------------
c
c..set matricies for eigenvalue equation
c
      do j1=1,idim
        zalfr(j1) = 0.0
        zalfi(j1) = 0.0
        zbeta(j1) = 0.0
        do j2=1,idim
          zamr(j1,j2) = 0.0
          zami(j1,j2) = 0.0
          zbmr(j1,j2) = 0.0
          zbmi(j1,j2) = 0.0
          zvr(j1,j2) = 0.0
          zvi(j1,j2) = 0.0
        enddo
      enddo
c
      if ( ieq .eq. 2 ) then
c
c..two equations when trapped particles and FLR effects omitted
c
        zamr(1,1) = ( 1.0 / zep2nh ) - ztauh - 1.0
        zamr(1,2) = - ztauh
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(1,1) = 1.0
        zbmr(1,2) = 0.0
        zbmr(2,1) = - tvr
        zbmr(2,2) = 1.0
c
      elseif ( ieq .eq. 4 ) then
c
c..4 equations with trapped electrons and FLR effects
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
      if  ( lprint .gt. 5 ) then
        write (6,*)
        write (6,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
      endif
c
c       ion continuity,with Singh RF terms 9/9 08.
c
        zamr(1,1) =1.0-zep2nh*rb2-zflh*rb3*ztauh*(1.0 + zetah)
        zami(1,1)=sigp
        zamr(1,2) = - zep2nh*rb2 * ztauh
        zami(1,2)=sigt
        zamr(1,3) = - zep2nh*rb2 * ztauh
c
        zbmr(1,1) = zflh*rb3*zep2nh
        zbmr(1,3) = zep2nh
c
c  ion energy
c
        zamr(2,1) = zetah - tvr
        zamr(2,2) = - zep2nh*rb2 * ztauh * ftr
c
        zbmr(2,2) =   zep2nh
        zbmr(2,3) = - zep2nh*tvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   by the ion density perturbation. The dilution factor 1-zfnz has now been 
c   added. Bounce averaged Omega_De included, zep2ne->zep2ne*(0.25+S*2/3)
c   !!!!!!!!Boumce average rbou included 25/10/05!!!!    
c   !!!  rbou=0.25d0+S*2.d0/3.d0
        zamr(3,1) = zft - zep2ne*rbou
        zamr(3,3) = zep2ne*rbou * (1. - zfnz - zfs)
        zamr(3,4) = zft * zep2ne*rbou
c
        zbmr(3,1) = ( zft - 1.0 ) * zep2ne
        zbmr(3,3) = zep2ne * (1. - zfnz -zfs)
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zetae - tvr )
        zamr(4,4) = zft * zep2ne*rbou * ftr
c
        zbmr(4,1) = ( 1.0 - zft ) * zep2ne * tvr
        zbmr(4,3) = - zep2ne * tvr
        zbmr(4,4) = zft * zep2ne
c
      else if ( ieq .eq. 6 ) then
c
c..Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
      if ( lprint .gt. 5 ) then
      write (6,*)
      write (6,*)
     & 'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      endif
c
c  hydrogen density
c
        zamr(1,1) = - 1.0
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.0
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.0
        zbmr(2,3) = - tvr
c
c  trapped electron density
c
        zamr(3,1) = - 1.0*rbou + zft/zep2ne
        zamr(3,3) = (1.0 - zfnz -zfs)*rbou
        zamr(3,4) = zft*rbou
        zamr(3,5) = zfnz*rbou
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz -zfs
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
c       zamr(4,1) och zamr(3,5) korrigerade 10/9/06
        zamr(4,1) = zft * ( zetae - tvr )/zep2ne
        zamr(4,4) = zft * ftr*rbou
c
        zbmr(4,1) = ( 1.0 - zft ) * tvr
        zbmr(4,3) = - ( 1.0 - zfnz -zfs) * tvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * tvr
c
c  impurity density
c
        zamr(5,1) = - 1.0
     &    + ( 1.0 - zflz * ztauz * ( 1.0 + zetaz ) ) / zep2nz
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = 1.0
c
c  impurity energy
c
        zamr(6,1) = ( zetaz - tvr ) / zep2nz
        zamr(6,6) = - ztauz * ftr
c
        zbmr(6,5) = - tvr
        zbmr(6,6) = 1.0
c
      else if ( ieq .eq. 7 ) then
c
c..Seven equations with impurities, trapped electrons, parallel ion motion
c and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here F is defined as F = GM*e phi/T_e where GM=1+etae/(epsn*(omega-1+i*vef))
c
      if ( lprint .gt. 5 ) then
      write (6,*)
      write (6,*)
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and Vp'
      endif
c
c********
c      H=0.5*DABS(S)/q
      H=0.
c********
c
c  hydrogen density
c
      zamr(1,1) = -1.*rb2
     & + (1. - zflh*rb3*ztauh*(1. + zetah ))/zep2nh 
      zami(1,1) = -H
      zamr(1,2) = - ztauh*rb2
      zami(1,2) = -ztauh*H
      zamr(1,3) = - ztauh*rb2+W0I
      zami(1,3) = -ztauh*H2
C
      zbmr(1,1) = zflh*rb3
      zbmr(1,3) = 1.
c
c hydrogen energy
c
      zamr(2,1) = (zetah - tvr )/zep2nh
      zamr(2,2) = - ztauh*rb2*ftr
c
      zbmr(2,2) = 1.
      zbmr(2,3) = -tvr
c
c trapped electron density
c
      zamr(3,1) = -1.*rbou + zft/(zep2ne)
      zami(3,1) = vef*(1. - zft)
      zamr(3,3) = (1. - zfnz -zfs)*rbou
      zami(3,3) = -vef*(1. -zfnz - zfs)
      zamr(3,4) = zft*rbou
      zamr(3,5) = zfnz*rbou
      zami(3,5) = -vef*zfnz
      zami(3,7) = vef*zft
c
      zbmr(3,1) = zft - 1.
      zbmr(3,3) = 1. - zfnz -zfs
      zbmr(3,5) = zfnz
c
c trapped electron energy
c
      zamr(4,1) = zft * (zetae - tvr )/(zep2ne)
      zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
      zami(4,3) = -vef*tvr*bt1*(1. -zfnz - zfs)
      zamr(4,4) = zft * ftr*rbou
      zami(4,5) = -vef*tvr*bt1*zfnz
      zami(4,7) = -ftr*vef*zft
c
      zbmr(4,1) = (1. - zft) * tvr
      zbmr(4,3) = - (1. - zfnz -zfs) * tvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * tvr
c
c impurity density
c
      zamr(5,1) = -1.*rb2
     & + ( 1. - zflz*rb3 * ztauz * (1. + zetaz ))/zep2nz
      zami(5,1)=-zimp*H/A
      zamr(5,5) = - ztauz*rb2+W0Z
      zami(5,5) = -zimp*ztauz*H/A
      zamr(5,6) = - ztauz*rb2
      zami(5,6) = -zimp*ztauz*H/A
c
      zbmr(5,1) = zflz*rb3
      zbmr(5,5) = 1.
c
c impurity energy
c
      zamr(6,1) = (zetaz - tvr)/zep2nz
      zamr(6,6) = - ztauz*rb2 * ftr
c
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
c
c  variable F
c
      zamr(7,1) = zetae/zep2ne - 1.
      zami(7,1) = vef
      zamr(7,7) = 1.
      zami(7,7) = -vef
c
      zbmr(7,1) = -1.
      zbmr(7,7) = 1.
c
      else if ( ieq .eq. 8 ) then
c
c..Eight equations with impurities, trapped electrons, parallel ion
c  motion collisions and FLR 
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Vp 
c  Here F = omega*gamma (collisions).
c
c      OBS! Nya termer för RF inlagda i jon-ekv 3/12/07.
      k1 = 0.25*q*q*zflh*SQRT(DABS((1.+zetah)*ztauh/((1.-zft)*zep2nh)))
      k2 = q*q*zflh*zflh*(1.+zetah)*ztauh/((1.-zft)*zep2nh)
      alp = 0.5*(k1+SQRT(k1+S*S*k2))
      alf = alp/(2.D0*zflh*q*q*betae*(1.D0 - zft))
c      kps = 0.5*SQRT(alp/zflh)/q
      kps = 0.289/(q*SQRT(zflh))
c      kps=0.D0
c
c
c hydrogen density
c
        zamr(1,1) = - 1.0*rb2+W0I*zflh*rb3
     &   + ( 1.0 - zflh*rb3 * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zamr(1,2) = - ztauh*rb2
        zamr(1,3) = - ztauh*rb2+W0I
        zamr(1,8) = kps
c
        zbmr(1,1) = zflh*rb3
        zbmr(1,3) = 1.0
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh*rb2 * ftr+W0I
        zamr(2,3)=-tvr*W0I
c
        zbmr(2,2) = 1.0
        zbmr(2,3) = - tvr
c
c  total electron density expressed in ion and imp densities
c
        zamr(3,1) = - 1.0*rbou + zft / zep2ne 
        zami(3,1) = vef*(1.-zft)
        zamr(3,3) = (1.0 - zfnz -zfs)*rbou
        zami(3,3) = -vef*(1. - zfnz -zfs)
        zamr(3,4) = zft*rbou
        zamr(3,5) = zfnz*rbou
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz -zfs
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft*(zetae - tvr)/zep2ne 
        zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
        zami(4,3) = -vef*tvr*bt1*(1.-zfnz -zfs)
        zamr(4,4) = zft * ftr*rbou
        zami(4,5) = -vef*tvr*bt1*zfnz
        zami(4,7) = -ftr*vef*zft
c
        zbmr(4,1) = ( 1.0 - zft ) *tvr
        zbmr(4,3) = - ( 1.0 - zfnz -zfs ) *tvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * tvr
c
c  impurity density
c
        zamr(5,1) = - 1.0*rb2
     &    + ( 1.0 - zflz*rb3 * ztauz * ( 1.0 + zetaz ) ) / zep2nz
        zamr(5,5) = - ztauz*rb2+W0Z
        zamr(5,6) = - ztauz*rb2
c
        zbmr(5,1) = zflz*rb3
        zbmr(5,5) = 1.0
c
c  impurity energy
c
        zamr(6,1) = ( zetaz - tvr ) / zep2nz
        zamr(6,6) = - ztauz*rb2 * ftr
c
        zbmr(6,5) = - tvr
        zbmr(6,6) = 1.0
c
c  variable F
c
         zamr(7,1) = zetae/zep2ne - 1. 
         zami(7,1) = vef
         zamr(7,7) = 1. 
         zami(7,7) = -vef
c
         zbmr(7,1) = -1.
         zbmr(7,7) = 1.
c
c     Parallel ion motion Vpi/Cs
c
         zamr(8,1) = kps
         zamr(8,2) = kps*ztauh
         zamr(8,3) = kps*ztauh
c
         zbmr(8,8) = 1.
c
      ELSE
      GO TO 08886
      ENDIF
      GOTO 08888
08886 CONTINUE
c
c--------------------------------------------------------------------------------
C      else if ( ieq .eq. 9 ) then
      IF(IEQ.NE.9) GO TO 08887
c
c..Nine  equations with impurities, trapped electrons, parallel ion
c  motion, collisions,  FLR , finite beta and parallel motion of impurities
c
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Av, K
c  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
c   magnetic vector potential.
c
c      WRITE(*,20061) WZJ(9),IK
20061 FORMAT(' Enter 9 Eq:s WZJ(9)=',2G11.3,' IK=',I13)

      EN=zep2nh
      ENI=1./EN
      EI=zetah
      EE=zetae
      TAU=1./ztauh
      TAUI=ztauh
      FL=zflh
      FT=zft
c
      GM=1./(1.-FT)
      BTA=FTRT*(GM+TAUI)
      H1=4.*q*q*FL
      XT=1./(1.+TAUI)
c
      ALA=2.*em*q*q*betae*2.*(1.+EE+TAUI*(1.+EI))/EN
      SH2=2.*S-1.+(KAPPA*(S-1.))**2
      SH=SQRT(SH2)
      H=0.5*ABS(SH)/q
      H2=IU*H
      GAV=1.
      ITS=0
      ITERA=1
c
c      ISEARCH=1
c      SEARCHMODE=2
c  If ISB=1 we use the strong ballooning approx (GAV=1)
c     IST is 1 at the first step. Then we use an analytical approx.
      IF(ISB.EQ.1) GO TO 20001
      IF(IST.NE.1) GO TO 800
c********************************************************************
c    Here the analytical approximation for WZ is calculated 
c    and stored in WZJ(IK)
c
      E1=FTRT*(1.+FL)-ENI+FL*TAUI*ENI*(1.+EI)+(GM+FTRT)*GAV
     &+H2*(1.+FTRT)
      E1=0.5*E1/(1.+FL)
      E=(TAUI*ENI*GM*(EI-TVR)+BTA)*(GAV+H2)-FTRT*ENI*(1.-FL*TAUI*(1.
     &+EI))
      E=E/(1.+FL)
      WZ1=-E1+SQRT(E1*E1-E)
      WZ2=-E1-SQRT(E1*E1-E)
      WZ=WZ1
      IF(IMAG(WZ2).GT.IMAG(WZ1)) WZ=WZ2
      WZI=IMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
C      WRITE(*,10021)EI,E,EN,ENI,GM,BTA,GAV,FTRT,H2
10021 FORMAT(' EI=',G11.4,' E=',2G11.4,' EN=',G11.4,' ENI=',G11.4,
     &' GM=',G11.4,' BTA=',G11.4,' GAV=',G11.4,' FTRT=',G11.4,
     &' H2=',2G11.4)
c
      WZJ(IK)=WZ
      WZH=EN*WZ
      IF(LPRINT.EQ.2)   WRITE(*,10001) WZH
10001 FORMAT(2X,' DISP9T WZ=',2G11.3)
c*********************************************************************
C--THE ITERATIVE PROCEEDURE FOR WZ STARTS HERE --
c
      IMET=0
  800 CONTINUE
      WZ=WZJ(IK)
c----------------------------------------------------
      IF(lprint.NE.2) GO TO 30006
      WRITE(*,00798) WZ,ITERA,ISEARCH,IM
00798 FORMAT(' WZ=',2G11.3,' ITERA=',I5,' ISEARCH=',I5,'   IM=',I5)
      IF(IMET.EQ.1) GO TO 30001
      IF(IMET.EQ.2) GO TO 30002
      WRITE(*,30003)
30003 FORMAT(' Average WZ')
      GO TO 30006
30001 WRITE(*,30004) 
30004 FORMAT(' Newton-Raps')
      GO TO 30006
30002 WRITE(*,30005)
30005 FORMAT(' Muller')
c
c--------------------------------------------------
30006 CONTINUE
c
      WZP=WZ
      ALPK=0.5*SH*SQRT(H1*XT*FL*(1.+TAUI*(1.+EI)/(EN*WZ)))
      IF(DREAL(ALPK).GE.0.) GOTO 801
      ALPK=-ALPK
  801 CONTINUE
      ALPC=-IU*ALPK
      ALPHA=-IU*ABS(SH)*q*FL
      XH=ABS(ALPHA/ALPC)
      ALPC=XH*ALPC
      R=2.*ABS(DREAL(WZ*ALPC))
      IF(R.LT.0.001) R=0.001    !! NEW 01.03.8
      HQ=2.*ALPC/H1
  802 GAV=(1.+0.5*S/R)*EXP(-0.25/R)
      GAV=GAV-0.5*ALA*(1.-EXP(-1./R))
      IF(GAV.LT.0.01) GAV=0.01
c
20001 CONTINUE
      IF(ISB.NE.1) GO TO 20002
      k1=0.25*q*q*zflh*DSQRT(DABS((1.+zetah)*ztauh/((1.-zft)*zep2nh)))
      k2=q*q*zflh*zflh*(1.+zetah)*ztauh/((1.-zft)*zep2nh)
      R=k1+SQRT(ABS(k1+S*S*k2))
      H=0.5*DABS(SH)/q
      HQ=-IU*H
      ITC=2
20002 CONTINUE
      alp=0.5*R
      IF(alp.LT.0.1) alp=0.1
      alf = alp/(2.D0*zflh*q*q*betae*(1.D0 - zft))
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
      RAV=1.+0.25*SH2/alp
c      WRITE(*,10002) alp,SH2,RAV
10002 FORMAT(2X,'alp=',G11.3,' SH2=',G11.3,' RAV=',G11.3)
c      WRITE(*,10003) XH,GAV,alf
10003 FORMAT(2X,'XH=',G11.3,' GAV=',G11.3,' alf=',G11.3)
c
c
c  *********
      HQR=DREAL(HQ)
      HQI=DIMAG(HQ)
c  *********
c--- WE NOW DEFINE THE MATRIX ELEMENTS --------------
c hydrogen density
c
        zamr(1,1) = - GAV*rb2+HQR
     &   + ( 1.0 - zflh*rb3 * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zami(1,1) = HQI
        zamr(1,2) = (HQR-GAV)*rb2*ztauh
        zami(1,2) = ztauh*HQI
        zamr(1,3) = (HQR-GAV)*rb2*ztauh
        zami(1,3) = ztauh*HQI
        zamr(1,8) = -em*ztauh*HQR*(1.+zetah)/(kpc*zep2nh)
        zami(1,8) = -em*ztauh*HQI*(1.+zetah)/(kpc*zep2nh)
        zamr(1,9) = -em*HQR/kpc
        zami(1,9) = -em*HQI/kpc
c
        zbmr(1,1) = zflh*rb3
        zbmr(1,3) = 1.
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr )/zep2nh
        zamr(2,2) = - ztauh*rb2 * ftr
c
        zbmr(2,2) = 1.
        zbmr(2,3) = - tvr
c
c  total electron density expressed in ion density and imp density
c
       zamr(3,1) = -1.*rbou + zft/(zep2ne)
       zami(3,1) = vef*(1.-zft)
       zamr(3,3) = (1. - zfnz -zfs)*rbou
       zami(3,3) = -vef*(1. - zfnz - zfs)
       zamr(3,4) = zft*rbou
       zamr(3,5) = zfnz*rbou
       zami(3,5) = -vef*zfnz
       zami(3,7) = vef*zft
       zamr(3,8) = -em*(1. - zft)/(kpc*zep2ne*rbou)
       zami(3,8) = em*(1.-zft)*vef/(kpc*zep2ne*rbou)
       zamr(3,9) = em*(1. - zft)*(1.+eni)/kpc
       zami(3,9) = -em*(1.-zft)*vef/kpc
c
      zbmr(3,1) = zft - 1.
      zbmr(3,3) = 1. - zfnz - zfs
      zbmr(3,5) = zfnz
      zbmr(3,9) = em*(1. - zft)/kpc
c
c  trapped electron energy
c
      zamr(4,1) = zft*(zetae - tvr)/(zep2ne)
      zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
      zami(4,3) = -vef*tvr*bt1*(1.-zfnz -zfs)
      zamr(4,4) = zft*ftr*rbou
      zami(4,5) = -vef*tvr*bt1*zfnz
      zami(4,7) = -ftr*vef*zft
c
      zbmr(4,1) = (1. - zft)*tvr
      zbmr(4,3) = -(1. - zfnz -zfs)*tvr
      zbmr(4,4) = zft
      zbmr(4,5) = -zfnz*tvr
c
c  impurity density
c
      zamr(5,1) = - GAV*rb2 +zimp*HQR/A
     & +(1. -zflz*rb3*ztauz*(1.+zetaz))/zep2nz
      zami(5,1) = zimp*HQI/A
      zamr(5,5) = (HQR*zimp/A-GAV)*rb2*ztauz
      zami(5,5) = zimp*ztauz*HQI/A
      zamr(5,6) = (HQR*zimp/A-GAV)*rb2*ztauz
      zami(5,6) = zimp*ztauz*HQI/A
      zamr(5,8) = -em*HQR*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
      zami(5,8) = -em*HQI*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
      zamr(5,9) = -em*HQR*zimp/(kpc*A)
      zami(5,9) = -em*HQI*zimp/(kpc*A)
c
      zbmr(5,1) = zflz*rb3
      zbmr(5,5) = 1.
c
c  impurity energy
c
      zamr(6,1) = (zetaz - tvr)/zep2nz
      zamr(6,6) = -ztauz*rb2*ftr
c
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
c
c  variable F
c
      zamr(7,1) = zetae/zep2ne - 1.
      zami(7,1) = vef
      zamr(7,7) = 1.
      zami(7,7) = -vef
c
      zbmr(7,1) = -1.
      zbmr(7,7) = 1.
c
c
c  electromagnetic parallel vectorpotential Av = e A_par/Te
c
      fft=(1.-zfnz-zfs)/(1.-zft)
      fzft=zfnz/(1.-zft)
      zamr(8,1) = em1*kpc*(eni+HQR*(fft+zimp*fzft/A))
      zami(8,1) = em1*HQI*(fft+zimp*fzft/A)*kpc
      zamr(8,2) = em1*HQR*ztauh*fft*kpc
      zami(8,2) = em1*HQI*ztauh*fft*kpc
      zamr(8,3) = em1*HQR*ztauh*fft*kpc
      zami(8,3) = em1*HQI*ztauh*fft*kpc
      zamr(8,5) = em1*HQR*zimp*ztauz*fzft*kpc/A
      zami(8,5) = em1*HQI*zimp*ztauz*fzft*kpc/A
      zamr(8,6) = em1*HQR*zimp*ztauz*fzft*kpc/A
      zami(8,6) = em1*HQI*zimp*ztauz*fzft*kpc/A
      zamr(8,8) = em1*((1.+zetae)*eni - alf*zflh*RAV)
     &-em1*HQR*(fft*ztauh*(1.+zetah)/zep2nh
     &+zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*A))
      zami(8,8) = -em1*HQI*(fft*ztauh*(1.+zetah)/zep2nh
     &+zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*A))
      zamr(8,9)= -em1*(eni+HQR*(fft+zimp*fzft/A))
      zami(8,9) = -em1*HQI*(fft+zimp*fzft/A)
c
      zbmr(8,1) = em1*kpc
      zbmr(8,8) = em1
      zbmr(8,9)= -em1
c
c     K = omega*Av
c
      zamr(9,9) = em1
c
      zbmr(9,8) = em1
c
      GO TO 08888
c      else
08887 CONTINUE
c
      write(6,*)
      write(6,*) ieq,' = ieq in sbrtn disp9t'
c      call abortb(6
c     &,'the value of ieq is wrong in sbrtn disp9t')
c
c       endif
c
08888 CONTINUE
c---------------------------------------------------------------
c  This part is common to all values of ieq
c..save copy of matrix which is over-written by sbrtn f02bjf
c
      do j2=1,ieq
        do j1=1,ieq
          zamrt(j1,j2) = zamr(j1,j2)
          zbmrt(j1,j2) = zbmr(j1,j2)
        enddo
      enddo
c
c..diagnostic output
c
      if ( lprint .gt. 9 ) then
        write (6,*)
        write (6,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (6,192) (zamr(j1,j2),j2=1,ieq)
        enddo
c
        write (6,*)
c
        write (6,*) ' zami(j1,j2)  j2 ->'
        do j1=1,ieq
          write(6,192) (zami(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
        write (6,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (6,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
       write(6,*) ' zbmi(j1,j2)  j2->'
       do j1=1,ieq
         write(6,192) (zbmi(j1,j2),j2=1,ieq) 
       enddo
c
 192  format (1p5e12.4)
      endif
c
c      write(*,193) zep2ne,zep2nh,zep2nz
  193 format(2X,'zep2n=',G12.5,' zep2nh=',G12.5,' zep2nz=',G12.5)
c      zgne=2./zep2ne
c      zgnh=2./zep2nh
c      zgnz=2./zep2nz
c
c      write(*,194) zgne,zgnh,zgnz
  194 format(2X,'zgne=',G12.5,' zgnh=',G12.5,' zgnz=',G12.5)
c

c..find the eigenvalues using NAG14 routine  f02gjf
c
      ztol = max ( 0.D0, cetain(32) )
      ifail = 1
c
  201 continue
c
      call r8tomsqz( idim,ieq,zamr,zami,zbmr,zbmi,
     &zalfr,zalfi,zbeta,zvr,zvi,ifail )
c
  202 continue
c
       if( ifail .le. 0 ) goto 210
       iret = 1
c      write(*,205) zetah,zetaz,zetae,ztauh,ztauz,zep2nh,zep2nz
c     &,zep2ne,zft,zfnz
c  205 format(x,'zetah',G10.4,' zetaz',G10.4,' zetae',G10.4,' ztauh',
c     1G10.4,/,' ztauz',G10.4,' zep2nh',G10.4,' zep2nz',G10.4,' zep2ne'
c     1,G10.4,/,' zft',G10.4,' zfnz',G10.4)
      go to 215
c
  210 continue
c..compute the complex eigenvalues
c
      do j=1,ieq
c
        zb = zbeta(j)
        if ( abs(zbeta(j)) .lt. zepsilon ) zb = zepsilon
        omega(j) = zalfr(j) / zb
        gamma(j) = zalfi(j) / zb
        ZZ(j)=omega(j)+(0.D0,1.D0)*gamma(j)
c
       enddo
c
  220 format(2X,'wr=',G11.4,' wi=',G11.4)
c
  215 continue
c
        if ( lprint .gt. 6 ) then
        do j=1,ieq
        write(*,220) omega(j),gamma(j)
      enddo
      endif
c--------------------------------------------------------------------
      IF(ieq.le.8) go to 00086 
      IF(ISB.EQ.1) go to 00086
      IF(ISEARCH.EQ.3) go to 10083
c**************************************************************************
c**  HERE THE PROCEEDURE TO FIND THE NEW WZ STARTS  ******************
c
c  Find the fastest growing mode, IM. If ISEARCH=1 only modes propagating
c  in the ion drift direction are considered.
      WM=0.001
      IM=0
      DO 00082 j=1,ieq
      IF(ISEARCH.EQ.2) GO TO 20021
      IF(omega(j).GE.0.) GO TO 00082
20021 IF(gamma(j).LT.WM) GO TO 00082
      WM=gamma(j)
      IM=j
00082 CONTINUE
C
      IMX=MAX(1,IM)
      WS=ZZ(IMX)
      IF(ITERA.EQ.1.AND.ISEARCH.EQ.1) WSF=WZP
c --  DEFINE THE NEW WZ ----
c -- If there is no unstable mode only average WZ is used--
      IMET=0
      IF(IM.GE.1) GO TO 10025
c      WZ=(WZ+WZP)/2.
      IF(ITERA.GT.1) GO TO 20019
      WZI=IMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
      WZA=WZ
      WZ=(WZ+WZP)/2.
      IF(lPRINT.EQ.2) WRITE(*,10028) WZA
10028 FORMAT(' WZA=',2G11.3)
      GO TO 00083
20019 CONTINUE
      WZ=(WZ+WZP)/2.
      GO TO 00083
10025 CONTINUE
      WZ=WS
c----------------------------------
c A new WZ corresponding to an unstable mode has been found --
c ----------------------------------
10026 CONTINUE
c  CHECK FOR CONVERGENCE 
c Independently of method used the average is always taken first
      WZ=(WZ+WZP)/2.
      IF(ITC.EQ.0) GO TO 00083
c      CTEST=ABS((WZ-WZP)/WZP)
      CTEST=ABS((WS-WSF)/WSF)
c      CTEST=ABS((WS-WZ)/WS)
      IF(CTEST.LE.TOL) ITS=1
      IF(ITS.EQ.1) GO TO 00083
c-----
      WZI=IMAG(WZ)
      IF(ITERA.GT.1) GO TO 10011
      WZI=IMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
      WZA=WZ
c -- In the first iteration we only use average
       GO TO 00083
c------------------------------------
10011 CONTINUE
c-- If ITC=2 we only use average ---------
      IF(ITC.EQ.2) GO TO 00083
c---------------------------------
c-- In the second and higher iteration we expand W=W(WK)
c-- In the second iteration we use Newton-Raphsons method--
      DWN=WZP-WZPP
      IF(ABS(DWN).LT.0.001) GO TO 00083
      DEW=(WZ-WZP)/DWN
c-- In the third and higher iterations we use Mullers method--
      IF(ITERA.GE.3) GO TO 10027  !! We will use Mullers method
c--------------------------------
20003 CONTINUE
      IF(ABS(1.-DEW).LT.0.001) GO TO 00083
      WZ=(WZ-WZP*DEW)/(1.-DEW)  !! This is the Newton Raphson result
      IMET=1
      GO TO 00083
10027 CONTINUE
c---- We will here use Mullers method -----------
      DWN1=WZ-WZP
      IF(ABS(DWN1).LT.0.001) GO TO 00083
c      DEW=(WSF-WSFF)/DWN
      D2EW=(DEW-DEWP)/DWN1
      IF(ABS(D2EW).LT.0.001) GO TO 20003
      AM1=(1.+D2EW*WZP-DEW)/D2EW
      BM1=2.*(WZP*DEW-WSF)/D2EW-WZP*WZP
      DR=SQRT(AM1**2+BM1)
      W1=AM1+DR
      W2=AM1-DR
      DW1=ABS(W1-WZP)
      DW2=ABS(W2-WZP)
      WZ=W1
      IF(DW1.GT.DW2) WZ=W2
c     We have now obtained WZ by Mullers method 
      IMET=2
c------------------------------------------------------------
00083 CONTINUE
c -- THE NEW WZ HAS BEEN OBTAINED
      WZI=IMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI) !! WZI always larger than 0.01
      WZI=IMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI) !! May  sometimes  be necessary  
      IF(ABS(WZ).GT.1000.) WZ=WZ*1000./ABS(WZ)  !! NEW 01.03.08
      WZJ(IK)=WZ
c      IF((ABS((WZ-WZP)/WZP)).LE.TOL) ITS=1
      IF(IM.EQ.0) ITS=0
c
C      WRITE(*,20032) IM,ITS,WZJ(IK),IK
20032 FORMAT(' IM=',I5,' ITS=',I5,' WZJ=',2G11.3,' IK=',I5)
c ---  Different conditions for continued iteration are tested ----
      IF(ITC.EQ.0) GO TO 00085
      IF(ITS.EQ.1) GO TO 00085
      IF(ITERA.GE.ITL) GO TO 00084
      IF(IM.EQ.0.AND.ITERA.GE.3) GO TO 00084
c------------------------------------------------------------
c  *** A new iteration will be made ***
20022 CONTINUE
      WZPP=WZP
      DEWP=DEW
      WZP=WZ
      WSFF=WSF
      WSF=WS
C
      DO j1=1,idim
      zalfr(j1)=0.
      zalfi(j1)=0.
      zbeta(j1)=0.
c
      DO j2=1,idim
      zamr(j1,j2)=0.
      zami(j1,j2)=0.
      zbmr(j1,j2)=0.
      zbmi(j1,j2)=0.
      zvr(j1,j2)=0.
      zvi(j1,j2)=0.
       enddo
      enddo
c
      ITERA=ITERA+1
      GO TO 800
c
00084 CONTINUE
      IF(IM.GE.1) GO TO 20024
      IF(ISEARCH.EQ.SEARCHMODE) GO TO 20024
      IF(ISEARCH.EQ.2) GO TO 20024
      ISEARCH=ISEARCH+1  !! Here the search is made with general WZ 
      WZJ(IK)=WZA
      ITERA=0
      GO TO 20022
20024 CONTINUE
C----WHEN THE ITERATIONS DO NOT CONVERGE WE USE THE FIRST AVERAGE--
      WZ=WZA
      WZR=DREAL(WZ)
      WZI=DIMAG(WZ)
      IF(WZR.LT.-10.) WZ=WZ-WZR-10.
      IF(WZI.GT.10.) WZ=WZ+IU*(10.-WZI)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
      WZI=IMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
C
      WZJ(IK)=WZ
      IF(lprint.NE.2) GO TO 00085
      WRITE(*,50001) WZ
50001 FORMAT(' NONCONVERGENT WZ=',2G11.4)
c----------------------------------------------------------------------
00085 CONTINUE
c New part using electron WZ for electron root 
      IF(SEARCHMODE.NE.3) GO TO 10085
C
      WRITE(*,90075)
90075 FORMAT(' SEARCMODE=3')
      IF(ISEARCH.EQ.3) GO TO 10083
      IF(DREAL(WZ).GT.0.) GO TO 10085
      ISEARCH=3
      IM=0
      WM=0.
      DO 10082 j=1,ieq
      WR=DREAL(ZZ(j))
      WI=DIMAG(ZZ(j))
      ZZS(j)=ZZ(j)
      IF(WR.LT.0.) GO TO 10082
      IF(WI.LE.WM) GO TO 10082
      IM=j
      WM=WI
10082 CONTINUE
c
      IF(WM.LT.0.001) GO TO 10085
      WZO=WZ
      WZ=ZZ(IM)
c
      DO j1=1,idim
      zalfr(j1)=0.
      zalfi(j1)=0.
      zbeta(j1)=0.
c
      DO j2=1,idim
      zamr(j1,j2)=0.
      zami(j1,j2)=0.
      zbmr(j1,j2)=0.
      zbmi(j1,j2)=0.
      zvr(j1,j2)=0.
      zvi(j1,j2)=0.
       enddo
      enddo
c
c      WRITE(*,10075) ISEARCH, WZ
10075 FORMAT(' ELECTRON MODE ISEARCH=',I5,' WZ=',2G11.3)
c
      GO TO 30006
10083 CONTINUE  !!! Here ISEARCH=3 and a new results has been obtained
c
      DO 10084 j=1,ieq
      IF(DREAL(ZZ(j)).GT.0.) GO TO 10084
      ZZ(j)=ZZS(j)
10084 CONTINUE
      WZE=WZ
      WZ=WZO
C----------------------------------------------------
10085 CONTINUE
      IF(lprint.NE.2) GO TO 00086
      WRITE(*,50002) EN*WZ,WZJ(IK),IK
50002 FORMAT(' EXIT DISP9T  WZ=',2G11.3,' WZJ(IK)=',2G11.3,'  IK=',I5)
      IF(SEARCHMODE.NE.3) GO TO 00086
      WRITE(*,10087) EN*WZE
10087 FORMAT(' SEARCHMODE=3, Electron WZ=',2G11.3)
00086 CONTINUE
      return
      end
