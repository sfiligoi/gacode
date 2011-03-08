c
        subroutine kapisn(
     >          nkimod,    !kapai model no. desired
     >          aimp,      !atomic mass of av. impurity
     >          xzimp,     !atomic number of av. impurity
     >          aplasm,    !array of atomic masses of hyd. species
     >          ng,        !number of hyd. species
     >          nj,        !max array size
     >          rhoel,     !zone centered electron density (cm**-3)
     >          rhob,      !z.c. array of hyd. spec. den_s (cm**-3)
     >          rhi,       !z.c. array of av. impurity_s density 
     >          rhoi,      !z.c. array of total ion density      
     >          te,        !z.c. array of Te (ev)
     >          ti,        !z.c. array of Ti (ev)
     >          zeff,      !z.c. array of plasma zeff
     >          numzones,  !number of zones 
     >          q,         !z.c. array of safety factor
     >          btf,       !tf (tesla) fld at cntr of outer flux 
     >          drshaf,    !shaf shift of outer zone bdr. (cm)
     >          rminor,    !plasma minor radius (cm)
     >          rmajor,    !major radius (cntr of outer flux) (cm)
     >          istringer, !stringer correction
     >          xkapi,     !z.c.  Neo-Class Ion thermal diff. (cm**2/sec) 
     >          xnstari,   !nu-star-ions, 1 at  banana-plateau transition
     >          xnstare,   !nu-star-elecs,1 at  banana-plateau transition 
     >          ztaui,     !ion collision time
     >          zrhoi,     !ion poloidal gyro-radius
     >          zfluxlim)  !flux limited max temp gradient length (cm)
c
        implicit none
c        include 'glf.m'
c
        integer ng, nj, numzones, nkimod, istringer, ih, j
        real*8 btf, rminor, rmajor, drshaf
        real*8 zk20, za2, zb2, zc2, zr0pr,
     >         zlamd, aratio, aratii, bpol, zmass, znhy,
     >         ztauh, zrhoh, zeps, zhstar, zalpha, zcoef,
     >         zkhbp, zkhpc, zkips, zchrg, zchrg4, ztau
        real*8 zvth, zeps32, znuhat, znstar, zrho, zrho2,
     >         rminloc, zbanana, zstring, zkb0b2, zk2st,
     >         zfps, zepsc, zb2b0i, zf, zhp, zepsq, zk2ht,
     >         zmstar, zk2, znhatb, zkcoef, zkbolt
        real*8 aplasm(ng), rhoel(nj), rhob(nj,3),
     >         rhoi(nj), te(nj), ti(nj), q(nj),
     >         xkapi(nj), xnstari(nj),
     >         xnstare(nj),rhi(nj),
     >         zeff(nj), aimp(nj), xzimp(nj),
     >         zfluxlim(nj), zrhoi(nj), ztaui(nj)

c       Dimension local arrays
c
        real*8 clogi(301)
c
c
c       2.67    CALCULATE NEO-CLASSICAL ION HEAT CONDUCTIVITY
c                       USING RUTHERFORD_S JULICH PAPER.
c               TO ACCOUNT FOR SEVERAL HYDROGENIC SPECIES, THE
c                       ATOMIC MASS HAS BEEN REDEFINED FOLLOWING
c                       RUTHERFORD_S APPROACH AS USED IN THE 
c                       BALDUR CODE
c
c  REVISED:  JUN-20-79  D. MC CUNE   NEOCLASSICAL CONDUCTIVITY
c          IS MULTIPLIED BY NAMELIST-INPUT FUDGE FACTOR XKFAC
c          AFTER IT IS CALCULATED. (XKFAC DEFAULTS TO 1.0)
c
c  REVISED:  JAN-6-81  D. MC CUNE  *** NEW COMMON VARIABLE ***
c         "NKIMOD" TO SELECT ION KAPA MODEL:  NKIMOD=1==>
c         OLD MODEL (RUTHERFORD--JULICH); NKIMOD=2==> MODIFIED
c         HAZELTINE-HINTON, AS USED IN BALDUR CODE.  THE
c         HAZELTINE-HINTON FORMULA FIT TO NEOCLASSICAL ION
c         THERMAL CONDUCTIVITY INTEGRALS WAS INTENDED FOR A
c         SINGLE SPECIES PLASMA; HERE IT IS USED APPROXIMATELY
c         GENERALIZED FOR MULTI-SPECIES PLASMA, FOLLOWING THE
c         APPROACH IN BALDUR; VALID ONLY IF ZEFF IS CLOSE TO
c         1.0.  (NOTE THESE COMMENTS ALSO APPLY TO JULICH FIT;
c         FOR SINGLE SPECIES ZEFF=1 PLASMA THE HAZELTINE-HINTON
c         FIT IS SUPERIOR).  THE RANGE OF VALIDITY OF THE
c         APPROXIMATIONS USED HEREIN IS A TOPIC WORTHY OF 
c         INVESTIGATION.
c            NOTE:  DEUTERIUM-HYDROGEN AND TRITIUM-DEUTERIUM
c         MIXING MAY ALSO EFFECT VALIDITY OF HAZELTINE-HINTON
c         (AND JULICH) APPROXIMATIONS
c  10-JULY 1981  D. MC CUNE
c    ADDING NKIMOD=3 OPTION TO USE C.W. BOLTON NEOCLASSICAL 
c  CONDUCTIVITY MODEL: SEE MAY 1981 U. OF TEXAS DISSERTATION
c  "THE ION THERMAL CONDUCTIVITY FOR A PURE TOKAMAK PLASMA"
c  BY C.W. BOLTON FOR HIS PHD.
c  ACCESSED BY SETTING INPUT VARIABLE NKIMOD=3.
c  NOTE BALDUR STYLE ZEFF APPROXIMATE- MODIFICATION IN COLLISIONALITY
c  IS STILL APPLIED AS AN ATTEMPT TO DEAL WITH IMPURE PLASMA, BUT 
c  THIS IS APPARENTLY INVALID IN THE CONTEXT OF THE BOLTON 
c  CALCULATION, THEREFORE NKIMOD=3 SHOULD BE USED ONLY IF ZEFF IS
c  SMALL *******
c
c  REV 10 SEP 1981  D. MC CUNE:
c   NEW PARAMETER XKFAC2 ADDED.  THIS ANOMOLOUS MULTIPLIER IS APPLIED
c  TO KAPA(I) AT EACH SURFACE WHERE Q .LE. 1.  
c
c   stan kaye -- adding chang-hinton kapa(i) (NKIMOD=4)
c     hazeltine hinton with shaf shift and finite epsilon corrections
c   PHYS. FLUIDS VOL. 25 SEPT 1982
c
c   DMC - NKIMOD=5:  CHANG HINTON WITH ADDED CORRECTIONS FOR ZEFF.NE.1
c    AS PER LETTER 20 JULY 1984 FROM DR CHANG TO DR GOLDSTON
c
c       KI*NI = (KH(BP) + KH(PS + CL)) * NH + KI(PS) * NIMP
c         (not sure what above comment refers to - RJG 9/7/87)
c
c   NKIMOD=5 model modified to agree with Chang and Hinton Physics
c       of Fluids article, (1986) p 3314
c       RJG 9/7/87
c
c-----------------------------------------------------------------
c       Modifications to make this a viable subroutine for SNAP
c       Eliminated anomalous factors in this code (all external now)
c       Note TRANSP version defined aratio as r/R, this one is
c       corrected in notation, but gives same physics
c       Rob Goldston 8/25/84
c
c       Modified to include Zeff(r).
c       Mike Zarnstorff   4/29/85
c-----------------------------------------------------------------
c
c  HAZELTINE-HINTON COEFFICIENTS
        data zk20/.66D0/
        data za2/1.03D0/
        data zb2/0.31D0/
        data zc2/0.74D0/
c
c-----------------------------------------------------------------
c
c       CALCULATE THE ION COULOMB LOGARITHM
c
        do 100 j=1,numzones
        zlamd=1.54D10*ti(j)/zeff(j)*dsqrt(te(j)/rhoel(j))
c       BOMB SHELTER: RJG 3/26/85
        if (zlamd .gt. 0.0) then
            clogi(j)=dlog(zlamd)
        else
            clogi(j)=17.D0
        endif
 100    continue
c
c       CALCULATE THE NEOCLASSICAL ION-HEAT CONDUCTIVITY
c
        do 200 j=1,numzones
c
        aratio = rmajor/(rminor*(j-0.5D0)/dfloat(numzones))
        aratii = 1.D0/aratio
        bpol = btf/(aratio*q(j))
C
C       AVERAGE HYDROGENIC MASS AND DENSITY
C
        zmass=0.D0
        znhy=0.D0
        do 210 ih=1,ng
          zmass=zmass+rhob(j,ih)*aplasm(ih)
          znhy=znhy+rhob(j,ih)
 210    continue
C
        zmass=zmass/znhy
C  BRANCH ACCORDING TO CHOICE OF CONDUCTIVITY MODEL
        go to (220,240,240,240,240) nkimod
 220    continue
C
C   JULICH  RUTHERFORD
C
        ztauh=1.473D7*dsqrt(zmass)*ti(j)**1.5D0/znhy/
     >        clogi(j)
        zrhoh=1.445D-2*dsqrt(zmass*ti(j))/(btf)
        zeps=(aratio)**1.5D0        !suspicious but as in TRANSP
        zhstar=4.897D-14*znhy*clogi(j)*aratii*
     >         btf*rmajor*zeps/(bpol*ti(j)**2.D0)
        zalpha=xzimp(j)*xzimp(j)*rhi(j)/(znhy)
        zcoef=zrhoh*zrhoh*q(j)**2.D0/ztauh
        zkhbp=zcoef*zeps*.48D0*(1.D0+1.5D0*zalpha)/
     >       (1.D0+.366D0*(1.D0+1.5D0*zalpha)*zhstar)
        zkhpc=zcoef*(.7D0/q(j)/q(j)+1.13D0+zalpha*
     >        (.5D0+.55D0/(.59D0+zalpha)))
        zkips=zcoef*dsqrt(aimp(j)/zmass)*rhi(j)/znhy*zalpha*1.13D0
c
        xkapi(j)=((zkhbp+zkhpc)*znhy+zkips*rhi(j))/rhoi(j)
c
c       Calculate nustari as in Hazeltine Hinton anyway

        zchrg=1.D0
        zchrg4=zchrg**4.D0
        ztau=2.083D7*dsqrt(zmass*ti(j)**3.D0)/
     >        (zchrg4*clogi(j))
c  IMPURITIES CORRECTION  (NORMAL HAZELTINE-HINTON WOULD HAVE
c  HYDROGENIC DENSITY IN DENOMINATOR)
        ztau=ztau/(zeff(j)*rhoel(j))
c
        zvth=dsqrt(3.204D-12*ti(j)/(1.67D-24*zmass))
        zeps=aratii
        zeps32=dsqrt(zeps**3.D0)
c
c  NU-STAR
c
        znuhat=1.414214D0*aratii*btf*rmajor/(bpol*zvth*ztau)
        znstar=znuhat/zeps32
        xnstari(j)=znstar
        xnstare(j)=znstar*1.414214D0*(ti(j)/te(j)/zchrg)**2.D0
c
        go to 200
c
c  HAZELTINE HINTON / BOLTON / chang hinton
c
 240    continue
c
        zchrg=1.D0
        zchrg4=zchrg**4.D0
        ztau=2.083D7*dsqrt(zmass*ti(j)**3.D0)/
     >        (zchrg4*clogi(j))
c  IMPURITIES CORRECTION  (NORMAL HAZELTINE-HINTON WOULD HAVE
c  HYDROGENIC DENSITY IN DENOMINATOR)
        ztau=ztau/(zeff(j)*rhoel(j))
c
        zvth=dsqrt(3.204D-12*ti(j)/(1.67D-24*zmass))
        zeps=aratii
        zeps32=dsqrt(zeps**3.D0)
c
c  NU-STAR
        znuhat=1.414214D0*aratii*btf*rmajor/(bpol*zvth*ztau)
        znstar=znuhat/zeps32
        xnstari(j)=znstar
        xnstare(j)=znstar*1.414214D0*(ti(j)/te(j)/zchrg)**2
c
        zrho=1.445D-2*dsqrt(zmass*ti(j))/(zchrg*bpol)
        zrho2=zrho**2.D0
        zrhoi(j)=zrho
        ztaui(j)=ztau

crew begin ********************************************
c 8.9.95 fix for  Stringer and flux limited flow
c 
c istringer=0 recovers old formula
c istringer=1 accounts for reduction of transport inside
c a banana width of r=0
c zfluxlim(j) is an output giving maximum allowed 
c temperature scale length -d ln Ti /d rho < 1/zfluxlim
c
        rminloc=rmajor*aratii+1.D-10
        zbanana=dsqrt(rminloc/rmajor)*dabs(zrho)
        zstring=zbanana**.66666D0*(rminloc)**.333D0
        if (istringer.eq.0) then
          zeps=rminloc/rmajor
        endif
        if (istringer.eq.1) then
          zeps=rminloc/rmajor
          if(rminloc.le.zbanana) zeps=zstring/rmajor
        endif
c
c now repeat calc for 
c
        zeps32=dsqrt(zeps**3.D0)
c
c  NU-STAR
c
        znuhat=1.414214D0*aratii*btf*rmajor/(bpol*zvth*ztau)
        znstar=znuhat/zeps32
        xnstari(j)=znstar
        xnstare(j)=znstar*1.414214D0*(ti(j)/te(j)/zchrg)**2
c
crew end ***********************************************
c
c  HAZELTINE-HINTON / chang-hinton:
        if (nkimod.eq.4) then
c  CHANG HINTON ORIGINAL VSN:
c
        if(j.eq.1) then
          zr0pr=0.D0
        else
c          zr0pr=(drshaf(j)-drshaf(j-1))/(rminor/dfloat(numzones))
          zr0pr=drshaf/(rminor/dfloat(numzones))
        endif
c
         zkb0b2=(1.D0+1.5D0*(zeps**2+zeps*zr0pr)+3.D0*zeps**3*
     >          zr0pr/8.D0)/(1.D0+0.5D0*zeps*zr0pr)
         zk2st=(1.D0+2.85D0*dsqrt(zeps)-2.33D0*zeps)*zkb0b2
         zfps=zeps32
        ENDIF
c
        IF(NKIMOD.EQ.5) THEN
c  CHANG HINTON WITH IMPROVED IMPURITY CORRECTION AS PER LETTER TO
c  R GOLDSTON 20 JULY 1984 FROM DR C.S. CHANG
c  THE ONLY IMPURITY CORRECTION IN THE ORIGINAL VSN WAS THE ENHANCEMENT
c  TO COLLISIONALITY AND ZTAU=THE ION-ION COLLISION TIME
c
c  Modified to follow Physics of Fluids (1986) p 3314 article by Chang
c       and Hinton, rather than private communication.
c
         if(j.eq.1) then
           zr0pr=0.D0
         else
c           zr0pr=(drshaf(j)-drshaf(j-1))/(rminor/dfloat(numzones))
           zr0pr=drshaf/(rminor/dfloat(numzones))
         endif
c
         zkb0b2=(1.D0+1.5D0*(zeps**2+zeps*zr0pr)+3.D0*zeps**3*
     >          zr0pr/8.D0)/(1.D0+0.5D0*zeps*zr0pr)
         zepsc=dsqrt(1.D0-zeps**2.D0)
         zb2b0i=zepsc*(1.D0+0.5D0*zeps*zr0pr)/
     >         (1.D0+zr0pr/zeps*(zepsc-1.D0))
         zf=0.5D0*(zkb0b2-zb2b0i)/dsqrt(zeps)
c        ZEFFAC=ZEFF(j)*RHOEL(J)/RHOI(J)
c        Looks like letter had a different impurity parameter
         zalpha = xzimp(j)**2*rhi(j)/znhy      !as for Rutherford/julich
c  PFIRSCH SCHLUTTER CORRECTION ("HP*F" IN THE LETTER):
c        ZHP=1.0+(1.0-1.0/zalpha)*(1.03/(zalpha-0.87)-0.68)
c        New version of Hp:
         zhp = 1.D0 + 1.33D0*zalpha*(1.D0+0.6D0*zalpha)/
     >         (1.D0+1.79D0*zalpha)
         zfps=zhp*zf
c  BANANA CORRECTION ("HB/0.66" IN THE LETTER)-- AGREES WITH ORIGINAL 
c  VSN IN LIMIT ZEFFAC=1.0
c        zk2st=(0.66+(1.88*dsqrt(zeps)-1.54*zeps)*(0.5D09+0.41/zalpha))
c     >         *zkb0b2/0.66
c        The guy above needs upgrading:
         zepsq = dsqrt(aratii)
         zepsq = dsqrt(zeps)
         zk2ht = ( 0.66D0*(1.D0+1.54D0*zalpha)+
     >           (1.88D0*zepsq-1.54D0*zeps)
     >           *(1.D0+3.75D0*zalpha) ) * zkb0b2
c        Chang-Hinton 1986 introduces a mu*star, and leaves tau as
c        the bulk-bulk collision frequency
         ztau = ztau * (zeff(j)*rhoel(j)) / znhy
         zmstar = znstar * (1.D0+1.54D0*zalpha)*znhy / 
     >            (zeff(j)*rhoel(j))
         endif
c
c  ORIGINAL HAZELTINE HINTON-- NO ADDL CORRECTIVE FACTORS:
        if (nkimod.eq.2) then
          zk2st=1.D0
          zfps=zeps32
        endif
c
c  HH AND DESCENDENT MODELS PUT TOGETHER:
c
        if ((nkimod.eq.2).or.(nkimod.eq.4)) then
          zk2=zk20*(zk2st/(1.D0+za2*dsqrt(znstar)+zb2*znstar) + zfps*
     >        zeps32*(zc2*zc2/zb2)*znstar/(1.D0+zc2*zeps32*znstar) )
          go to 245
        endif
        if (nkimod.eq.5) then
          zk2=zk20*(zk2ht/0.66D0/(1.D0+za2*dsqrt(zmstar)+
     >        zb2*zmstar)+zfps*
     >    zeps32*(zc2*zc2/zb2)*zmstar/(1.D0+zc2*zeps32*zmstar) )
          go to 245
        endif
c
c  BOLTON
c
        znhatb=.75D0*dsqrt(3.1415926D0)*znuhat
c  BOLTON_S COLLISIONALITY DIFFERS FROM HAZELTINE-HINTON_S BY FACTOR
c  3*SQRT(PI)/4 ...
        zk2=zkbolt(znhatb,zeps)
c
 245    zkcoef=zrho2/ztau
        xkapi(j)=zk2*dsqrt(aratii)*zkcoef
c
crew *******************************************************
        xkapi(j)=zk2/ztau*
     >     dsqrt(zeps)*(zbanana+1.D-10)**2/(zeps+1.D-10)
        if (istringer.eq.1..and.rminloc.le.zbanana) then
          xkapi(j)=xkapi(j)*rminloc/(zbanana+1.D-10)
        endif
        zfluxlim(j)=zbanana+1.D-10
        if (rminloc.le.zbanana) zfluxlim(j)=zstring+1.0D-10
crew *******************************************************
c  
        go to 200
c
 200    continue
c
        return
        end
