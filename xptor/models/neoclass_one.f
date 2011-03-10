        subroutine kapisn_one(
     >          nkimod,    !model number desired
     >          aimp,      !atomic mass of avg. impurity
     >          xzimp,     !atomic number of avg. impurity
     >          aplasm,    !array of atomic masses of hyd. species
     >          ng,        !number of hydrogen species
     >          rhoel,     !zone centered electron density (cm**-3)
     >          rhob,      !hyd. spec. den_s (cm**-3)
     >          rhi,       !avg. impurity density 
     >          rhoi,      !total ion density      
     >          te,        !Te (ev)
     >          ti,        !Ti (ev)
     >          zeff,      !plasma zeff
     >          q,         !array of safety factor
     >          btf,       !tf (tesla) fld at cntr of outer flux 
     >          dshafdr,   !derivative of shaf shift of outer zone bdr.
     >          rminor,    !plasma minor radius (cm)
     >          rmajor,    !major radius (cntr of outer flux) (cm)
     >          istringer, !stringer correction
     >          xkapi,     !Neoclassical ion thermal diffusivity (cm**2/s)
     >          xnstari,   !nu-star-ions, 1 at  banana-plateau transition
     >          xnstare,   !nu-star-elecs,1 at  banana-plateau transition 
     >          ztau,      !ion collision time
     >          zrho,      !ion poloidal gyro-radius
     >          zfluxlim)  !flux limited max temp gradient length (cm)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  v2.67    CALCULATES NEO-CLASSICAL ION HEAT CONDUCTIVITY
c           TO ACCOUNT FOR SEVERAL HYDROGENIC SPECIES, THE
c           ATOMIC MASS HAS BEEN REDEFINED FOLLOWING
c           RUTHERFORD'S APPROACH USED IN THE BALDUR CODE
c
c... This is the one-gridpoint version of kapisn.
c... Adapted by J. Kinsey July 20, 1999
c
c... Declare variables
c
        implicit none
c
        integer ng, nkimod, istringer, ih
        real*8 aplasm(ng), rhob(ng,1),
     >    rhoel, rhoi, te, ti, q,
     >    dshafdr, xkapi, xnstari,
     >    xnstare, rhi, zeff, aimp, xzimp,
     >    zfluxlim, clogi, btf, rminor, rmajor
        real*8 zpi, zpmass, zone, zhalf, z3half, zthird, z2third
c
        real*8 zk20, za2, zb2, zc2,
     >    zlamd, aratio, aratii, bpol, zmass, znhy,
     >    ztauh, zrhoh, zeps, zhstar, zalpha, zcoef,
     >    zkhbp, zkhpc, zkips, zchrg, zchrg4, ztau,
     >    zvth, zeps32, znuhat, znstar, zrho, zrho2,
     >    rminloc, zbanana, zstring, zkb0b2, zk2st,
     >    zfps, zepsc, zb2b0i, zf, zhp, zepsq, zk2ht,
     >    zmstar, zk2, znhatb, zkcoef, zkbolt
c
c... Revision History
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
c
c  HAZELTINE-HINTON COEFFICIENTS
c
        data zk20/.66D0/
        data za2/1.03D0/
        data zb2/0.31D0/
        data zc2/0.74D0/
c
c-----------------------------------------------------------------
c
        zpi=dabs(dacos(-1.D0))
        zpmass=1.67D-24
        zone=1.D0
        zthird=1.D0/3.D0
        z2third=2.D0/3.D0
        zhalf=1.D0/2.D0
        z3half=3.D0/2.D0
c
c       CALCULATE THE ION COULOMB LOGARITHM
c
        zlamd=1.54D10*ti/zeff*dsqrt(te/rhoel)
c       BOMB SHELTER: RJG 3/26/85
        if (zlamd .gt. 0.0) then
            clogi=dlog(zlamd)
        else
            clogi=17.D0
        endif
c
c       CALCULATE THE NEOCLASSICAL ION-HEAT CONDUCTIVITY
c
        aratio = rmajor/dmax1(1.D-10,rminor)
        aratii = zone/aratio
        bpol = btf/(aratio*q)
c
c       AVERAGE HYDROGENIC MASS AND DENSITY
c
        zmass=0.D0
        znhy=0.D0
        do 210 ih=1,ng
          zmass=zmass+rhob(ih,1)*aplasm(ih)
          znhy=znhy+rhob(ih,1)
 210    continue
c
        zmass=zmass/znhy
c
c-----------------------------------------------------------------
c...Rutherford model
c
        if (nkimod.eq.1) then
c
        ztauh=1.473D7*dsqrt(zmass)*ti**z3half/znhy/clogi
        zrhoh=1.445D-2*dsqrt(zmass*ti)/btf
        zeps=(aratio)**z3half        !suspicious but as in TRANSP
        zhstar=4.897D-14*znhy*clogi*aratii*
     >         btf*rmajor*zeps/(bpol*ti**2.D0)
        zalpha=xzimp*xzimp*rhi/znhy
        zcoef=zrhoh*zrhoh*q**2.D0/ztauh
        zkhbp=zcoef*zeps*.48D0*(zone+z3half*zalpha)/
     >       (zone+.366D0*(zone+z3half*zalpha)*zhstar)
        zkhpc=zcoef*(.7D0/q/q+1.13D0+zalpha*
     >        (.5D0+.55D0/(.59D0+zalpha)))
        zkips=zcoef*dsqrt(aimp/zmass)*rhi/znhy*zalpha*1.13D0
c
        xkapi=((zkhbp+zkhpc)*znhy+zkips*rhi)/rhoi
c
c       Calculate nustari as in Hazeltine Hinton anyway

        zchrg=zone
        zchrg4=zchrg**4.D0
        ztau=2.083D7*dsqrt(zmass*ti**3.D0)/
     >        (zchrg4*clogi)
c  IMPURITIES CORRECTION  (NORMAL HAZELTINE-HINTON WOULD HAVE
c  HYDROGENIC DENSITY IN DENOMINATOR)
        ztau=ztau/(zeff*rhoel)
        zvth=dsqrt(3.204D-12*ti/(zpmass*zmass))
        zeps=aratii
        zeps32=dsqrt(zeps**3.D0)
c
c  NU-STAR
c
        znuhat=dsqrt(2.D0)*aratii*btf*rmajor/(bpol*zvth*ztau)
        znstar=znuhat/zeps32
        xnstari=znstar
        xnstare=znstar*dsqrt(2.D0)*(ti/te/zchrg)**2.D0
c
        endif
c
c-----------------------------------------------------------------
c...Hinton-Hazeltine / Bolton / Chang-Hinton models
c
        if (nkimod.eq.2 .or. nkimod.eq.3 .or. nkimod.eq.4 
     >     .or. nkimod.eq.5) then
c
        zchrg=zone
        zchrg4=zchrg**4.D0
        ztau=2.083D7*dsqrt(zmass*ti**3.D0)/
     >        (zchrg4*clogi)
c  IMPURITIES CORRECTION  (NORMAL HAZELTINE-HINTON WOULD HAVE
c  HYDROGENIC DENSITY IN DENOMINATOR)
        ztau=ztau/(zeff*rhoel)
        zvth=dsqrt(3.204D-12*ti/(zpmass*zmass))
        zeps=aratii
        zeps32=dsqrt(zeps**3.D0)
c
c  NU-STAR
        znuhat=dsqrt(2.D0)*aratii*btf*rmajor/(bpol*zvth*ztau)
        znstar=znuhat/zeps32
        xnstari=znstar
        xnstare=znstar*dsqrt(2.D0)*(ti/te/zchrg)**2
c
        zrho=1.445D-2*dsqrt(zmass*ti)/(zchrg*bpol)
        zrho2=zrho**2.D0

c begin ********************************************
c 8.9.95 fix for  Stringer and flux limited flow
c 
c istringer=0 recovers old formula
c istringer=1 accounts for reduction of transport inside
c a banana width of r=0
c zfluxlim is an output giving maximum allowed 
c temperature scale length -d ln Ti /d rho < 1/zfluxlim

        rminloc=rmajor*aratii+1.D-10
        zbanana=dsqrt(rminloc/rmajor)*dabs(zrho)
        zstring=zbanana**z2third*(rminloc)**zthird
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
        znuhat=dsqrt(2.D0)*aratii*btf*rmajor/(bpol*zvth*ztau)
        znstar=znuhat/zeps32
        xnstari=znstar
        xnstare=znstar*dsqrt(2.D0)*(ti/te/zchrg)**2
c
c-----------------------------------------------------------------
c...Chang-Hinton models
c
        if (nkimod.eq.4) then
c  CHANG-HINTON ORIGINAL:
c
         zkb0b2=(zone+z3half*(zeps**2+zeps*dshafdr)+
     >          3.D0*zeps**3*dshafdr/8.D0)/
     >         (zone+zhalf*zeps*dshafdr)
         zk2st=(zone+2.85D0*dsqrt(zeps)-2.33D0*zeps)*zkb0b2
         zfps=zeps32
         endif
c
        if (nkimod.eq.5) then
c  CHANG HINTON WITH IMPROVED IMPURITY CORRECTION AS PER LETTER TO
c  R GOLDSTON 20 JULY 1984 FROM DR C.S. CHANG
c  THE ONLY IMPURITY CORRECTION IN THE ORIGINAL VSN WAS THE ENHANCEMENT
c  TO COLLISIONALITY AND ZTAU=THE ION-ION COLLISION TIME
c
c  Modified to follow Physics of Fluids (1986) p 3314 article by Chang
c       and Hinton, rather than private communication.
c
         zkb0b2=(zone+z3half*(zeps**2+zeps*dshafdr)+
     >          3.D0*zeps**3*dshafdr/8.D0)/
     >         (zone+zhalf*zeps*dshafdr)
         zepsc=dsqrt(zone-zeps**2.D0)
         zb2b0i=zepsc*(zone+zhalf*ZEPS*dshafdr)/
     >         (zone+dshafdr/zeps*(zepsc-zone))
         zf=zhalf*(zkb0b2-zb2b0i)/dsqrt(zeps)
c        ZEFFAC=ZEFF*RHOEL/RHOI
c        Looks like letter had a different impurity parameter
         zalpha = xzimp**2*rhi/znhy      !as for Rutherford/julich
c  PFIRSCH-SCHLUTER CORRECTION ("HP*F" IN THE LETTER):
c        ZHP=1.D0+(1.D0-1.D0/zalpha)*(1.03D0/(zalpha-0.87D0)-0.68D0)
c        New version of Hp:
         zhp = zone + 1.33D0*zalpha*(zone+0.6D0*zalpha)/
     >         (zone+1.79D0*zalpha)
         zfps=zhp*zf
c  BANANA CORRECTION ("HB/0.66" IN THE LETTER)-- AGREES WITH ORIGINAL 
c  VSN IN LIMIT ZEFFAC=1.0
c        zk2st=(0.66D0+(1.88D0*sqrt(zeps)-1.54D0*zeps)*
c              (0.59D0+0.41D0/zalpha))
c     >         *zkb0b2/0.66D0
c        The guy above needs upgrading:
         zepsq = dsqrt(aratii)
         zepsq = dsqrt(zeps)
         zk2ht = ( 0.66D0*(zone+1.54D0*zalpha)+
     >           (1.88D0*zepsq-1.54D0*zeps)
     >           *(zone+3.75D0*zalpha) ) * zkb0b2
c        Chang-Hinton 1986 introduces a mu*star, and leaves tau as
c        the bulk-bulk collision frequency
         ztau = ztau * (zeff*rhoel) / znhy
         zmstar = znstar * (zone+1.54D0*zalpha)*znhy / 
     >            (zeff*rhoel)
         endif
c
c-----------------------------------------------------------------
c  ORIGINAL HINTON-HAZELTINE MODEL -- NO ADDL CORRECTIVE FACTORS:
        if (nkimod.eq.2) then
          zk2st=zone
          zfps=zeps32
        endif
c
c  HH AND DESCENDENT MODELS PUT TOGETHER:
c
        if ((nkimod.eq.2).or.(nkimod.eq.4)) then
          zk2=zk20*(zk2st/(zone+za2*dsqrt(znstar)+zb2*znstar) + zfps*
     >        zeps32*(zc2*zc2/zb2)*znstar/(zone+zc2*zeps32*znstar) )
        endif
        if (nkimod.eq.5) then
          zk2=zk20*(zk2ht/0.66D0/(zone+za2*dsqrt(zmstar)+
     >        zb2*zmstar)+zfps*
     >    zeps32*(zc2*zc2/zb2)*zmstar/(zone+zc2*zeps32*zmstar) )
        endif
c
c-----------------------------------------------------------------
c...Bolton model
c   Note: BOLTON COLLISIONALITY DIFFERS FROM HAZELTINE-HINTON
c   BY FACTOR 3*SQRT(PI)/4 ...
        if (nkimod.eq.3) then
          znhatb=.75D0*dsqrt(zpi)*znuhat
          zk2=zkbolt(znhatb,zeps)
        endif
c-----------------------------------------------------------------
        zkcoef=zrho2/ztau
        xkapi=zk2*dsqrt(aratii)*zkcoef
c
c*******************************************************
        xkapi=zk2/ztau*
     >     sqrt(zeps)*(zbanana+1.D-10)**2/(zeps+1.D-10)
        if (istringer.eq.1..and.rminloc.le.zbanana) then
         xkapi=xkapi*rminloc/(zbanana+1.D-10)
        endif
        zfluxlim=zbanana+1.D-10
        if (rminloc.le.zbanana) zfluxlim=zstring+1.0D-10
c*******************************************************
c
        else
          write(6,*) 'Error: invalid option for nkimod'
          return
        endif
c
        return
c
        end
c
c*******************************************************
c
c  BOLTONS KI(NU-HAT,EPS) FUNCTION (THESIS PP 74--78)
c
        real*8 function zkbolt(znuhat,zeps)
c
        implicit none
c
        real*8 znuhat, zeps, zse, zeps32, ztanc,
     >    ztbnc, ztcnc, ztdnc, znh, zkinc,
     >    ztaps, ztbps, ztcps, ztdps, zkips
        real*8 zone
c
c  N-C COEFFICIENTS FROM TABLE VI
        real*8 zanc(3),zbnc(3),zcnc(3),zdnc(3)
c  P-S COEFFICIENTS FROM TABLE VI
        real*8 zaps(3),zbps(3),zcps(3),zdps(3)
c
        data zanc/2.441D0,-3.87D0,2.19D0/
        data zbnc/.9362D0,-3.109D0,4.087D0/
        data zcnc/.241D0,3.40D0,-2.54D0/
        data zdnc/.2664D0,-.352D0,.44D0/
c
        data zaps/.364D0,-2.76D0,2.21D0/
        data zbps/.553D0,2.41D0,-3.42D0/
        data zcps/1.18D0,.292D0,1.07D0/
        data zdps/.0188D0,.180D0,-.127D0/
c
        zone=1.D0
        zse=dsqrt(zeps)
        zeps32=zse**3.D0
        ztanc=.66D0+((zanc(3)*zse+zanc(2))*zse+zanc(1))*zse
        ztbnc=((zbnc(3)*zeps+zbnc(2))*zeps+zbnc(1))/(zeps**.75D0)
        ztcnc=((zcnc(3)*zeps+zcnc(2))*zeps+zcnc(1))/(zeps32)
        ztdnc=((zdnc(3)*zeps+zdnc(2))*zeps+zdnc(1))/(zeps32)
c
c  N-C PART
c
        znh=dsqrt(znuhat)
        zkinc=ztanc/
     >  (1.D0+ztbnc*znh+ztcnc*znuhat+ztdnc*znuhat*znuhat)
c
        ztaps=((zaps(3)*zse+zaps(2))*zse+zaps(1))*zeps32
        ztbps=((zbps(3)*zse+zbps(2))*zse+zbps(1))*zeps32
        ztcps=((zcps(3)*zse+zcps(2))*zse+zcps(1))
        ztdps=((zdps(3)*zse+zdps(2))*zse+zdps(1))
c
c  P-S PART
c
        zkips=1.57D0*zeps32 + (ztaps+ztbps*znh)/
     >        (zone+ztcps*znuhat**1.5D0+ztdps*znuhat**2.5D0)
c
        zkbolt=zkinc+zkips
c
      return
      end
