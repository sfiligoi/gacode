  MODULE thermalization

      USE nrtype,                          ONLY : DP,I4B

      USE constnts,                        ONLY : xmassp,xmasse,charg4, &
                                                  charge
                                                  
      USE common_constants,                ONLY : PI, ROOT2PI,kevperg,  &
                                                  zeroc,izero,Ln_2,     &
                                                  kevpjou

      USE error_handler,                   ONLY : terminate,lerrno,iomaxerr

      USE io_gcnmp,                        ONLY : nlog

      REAL(DP)  vcvo, tstcx, emzrat        ! replaces common gecom


    CONTAINS


      SUBROUTINE fast_ion_parms (atw,atwf,ene,en,enn,ezero,ibcur,ibcx,ifirst,kj, &
                        nj,nion,nneu,pzero,te,vz,zsq,zsqf,bke,bki,ecrit,         &
                        emzrat,encap,ge,gi,gth,taus,fionx,rtstcx)
! ----------------------------------------------------------------------
!     this is essentially slow1 in Onetwo 
!     This subroutine evaluates the transfer functions N, Ge, Gi, Ke, and
!     Ki defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974) to
!     describe fast ion slowing down.  Ge and Gi are the fractions of
!     the fast ion energy transferred to the electrons and thermal ions,
!     respectively.  The remaining energy fraction is assumed lost due
!     to secondary charge exchange, i.e., charge exchange between fast
!     ions and thermal neutrals.  Similarly, Ke and Ki are the fractions
!     of the fast ion parallel momentum transferred to electrons and
!     thermal ions.  Moreover, N, Ge, and Ke are proportional to the
!     slowing-down times for particles, energy, and parallel momentum
!     and are used in subroutine slow2.  The fast ions may be beam ions
!     from neutral beam heating or alpha particles from fusion.
!     This subroutine has been modified to account for plasma rotation.
!     Expressions contained here are valid for Vrot/Vthi< = 1 (although
!     assumptions of Maxwellian target distributions may be invalid).
!     An additional factor, Gth, describes momentum and energy sources
!     to the target ions due to nonzero kinetic energy at thermalization
!     with the rotating target.  (gts: 14jul89)
! EXTERNALS
!     cxrv_m
!     cxr_m
!     bkef_m
!     gef_m
!     encapf_m
! INPUT
!     atw(1:nion)  mass number of(all)  ion species in plasma
!     atwf   mass number of fast iosn
!     nj           rho grid size
!     nion         total number of ion species
!     nneu         number neutrals (1 or 2) If 2 we take average desnity
!     ifirst       if 1 assumes some calcualtions below are 
!                  allready set on first call ( this rotuien is called
!                  in a loop over beams and energies per beam), loop over
!                  spatial grid is done below.
!     kj           for 2d arrays need the declared dimension of firt index
!     pzero(*)     mass fast ion * initial speed fast ion
!     vz           toroidal rotation speed (approximated as omega *<R^2>/<R>
!     enn(1:nj,1:nneu)   neutral density (1/cm^3)
!     ene(1:nj)          electron density
!     en(1:nj,1:nion)    ion densities
!     rtstcx       set to 1.0 
!     fionx        fionx allows testing of non-classical slowing down
!     ezero        beam energy,kev (full half or third as appropriate)
!     te(1:nj)     electron temp in kev!
!     zsqf         charge sq of fast ion
!     xmasse       electron mass grams
!     xmassp       ion      mass grams
!
! OUTPUT (all through argument list)
!     ecrit(1:nj)  fast ion energy(kev) at which 50% of KE energy is
!                  delivered  to ions,and 50% to electrons
!     emzrat(1:nj) the term sum[( Z*n)/(atwf*sum(z^2*n/atw)] see Callen paper
!     taus(1:nj)   fast ion slowing down time
!     encap(1:nj)
!     ge(1:nj)
!     gi(1:nj)
!     gth(1:nj)
!     bke(1:nj)
!     bki(1:nj)
! ----------------------------------------------------------------------
!
      USE nrtype,                          ONLY : DP,I4B


      IMPLICIT NONE

      INTEGER(I4B) j,k,kj,nj,ifirst,ibcx,ibcur,ivbcxit,nion,nneu

      REAL(DP)  atw(*), ene(*), en(kj,*), enn(kj,*), pzero(*), te(*), &
                vz(*), zsq(kj,*)
      REAL(DP)  bke(*), bki(*), ecrit(*), emzrat(*), encap(*), ge(*), &
                gi(*), gth(*), taus(*)

      REAL(DP) atwf,econst,tconst,tiny,sum1,sum2,tegt0,teerg,xlam,    &
               ztaue,erot,prot,erl0,prel0,vrat,taurat,vbeamcx,        &
               erel0avg,taurat_part,rtstcx,vfast_avg,vcxerr,          &
               fionx,ft,prat,ezero,erel0,zsqf
               

!
!  calculate constants
!
      econst = 0.5_DP * (4.5_DP * pi * xmassp/xmasse) ** 0.333333_DP  ! nominaly
      tconst = atwf * xmassp / (zsqf * xmasse)
      tiny   = 0.0001_DP

   DO j=1,nj !  begin loop over mesh points
        IF (ifirst .EQ. 0)  go to 20
!
!       calculate ecrit, the "critical energy" at which fast ions
!       transfer equal energy to electrons and thermal ions;
!       also calculate emzrat = mi*<Z>/(mf*[Z])
!
        sum1 = zeroc
        sum2 = zeroc
        DO  k=1,nion
           sum1 = sum1 + en(j,k)*zsq(j,k)
           sum2 = sum2 + en(j,k)*zsq(j,k)/atw(k)
        ENDDO
        ecrit(j)  = econst*te(j)*atwf*(sum2/ene(j))**0.666667        !ke
        emzrat(j) = sum1 / (atwf*sum2)
!
!       calculate taus, the Spitzer momentum exchange time
!       for electron-ion collisions
!
        tegt0   = MAX (te(j), tiny)  ! keV, avoid probs near plasma edge
        teerg   = 1.6e-9*tegt0
        xlam    = 24.0 - LOG (SQRT (ene(j))/(1.0e3*tegt0))
        ztaue   = SQRT (xmasse*teerg**3) &
                / (1.33333*root2pi*ene(j)*charg4*xlam)
        taus(j) = tconst*ztaue

!
!       calculate erel0, the fast-ion initial energy in the rotating
!       frame; this is  0.5 * mass_fast_ion * (v_fastion-v_rot)**2
!       = 0.5*mf*(vf**2-2vf*vrot+ vrot**2)
!       = 0.5*mf*vf**2 -mf*vf*vrot +0.5*mt*vrot**2 
!       NOTE that it is assumed that mt ( mass thermal)= mf( mass,fast ) ion
!       = ezero -mf*vf*vrot + 0.5*mf*vrot**2
!       = ezero -erot, with erot = -0.5*mf*vrot**2 +mf*vf*vrot
   20   erot  = zeroc
        prot  = zeroc
        IF (vz(1) .NE. zeroc) THEN
          erot  = (pzero(j)*vz(j) - 0.5 * atwf*xmassp*vz(j)**2)*kevperg
          prot  = atwf*xmassp*vz(j)
        END IF
        erel0 = ABS (ezero - erot) 

  
        prel0 = pzero(j) - prot
!
!       calculate vrat = SQRT (ecrit/ezero) and taurat = taus/taucx,
!       where taucx is the charge exchange time
!
        vrat   = SQRT (ecrit(j)/erel0)


        taurat = zeroc
        IF (ibcx .EQ. 0  )  go to 50
        IF (atwf .GT. 3.0)  go to 50
!
        IF(nneu == 2)THEN
                  taurat_part = taus(j) * (enn(j,1)+enn(j,2))
        ELSE
                  taurat_part = taus(j) * enn(j,1)
        ENDIF
        ivbcxit     = 0
        vbeamcx     = SQRT (2.0*erel0 /(kevperg*atwf*xmassp))
        IF (rtstcx .GT. -5.0) &
          vbeamcx = ABS (rtstcx) * vbeamcx   ! use rtstcx for testing
        erel0avg  = erel0
 450    IF (rtstcx .GT. zeroc) THEN
          taurat  = taurat_part * cxr_m(erel0/atwf)
          taurat  = taurat * rtstcx ! rtstcx is input, defaulted to 1.0
        ELSE ! use sigma(v)*v rather than Maxwellian average
!              assume neutrals have same bulk speed as thermal ions
!              i.e., use erel0 as relative energy
          taurat  = taurat_part * cxrv_m(erel0avg)*vbeamcx
        END IF
!
!       calculate N( = encap), ge, and gi
!       assume taurat independent of v:
!

   50   encap(j) = encapf_m(vrat,taurat)

        ge(j)    = gef_m(vrat,taurat)*erel0/ezero
        IF (ge(j) .LT. zeroc)  ge(j) = zeroc
        IF(ge(j) .LE. zeroc)ge(j) = 1.e-15_DP !HSJ 4/6/01
!
        IF (rtstcx .LT. -5.0_DP) THEN ! use average fast ion speed
!                                    valid only after beam has
!                                    established itself .. HSJ .. 6/5/96
          erel0avg  = erel0*ge(j)/(2.*encap(j))
          vfast_avg = SQRT (erel0avg/(kevperg*atwf*xmassp))
          vcxerr    = ABS ((vfast_avg-vbeamcx)/vfast_avg)
          vbeamcx   = 0.5*(vfast_avg+vbeamcx)
          ivbcxit   = ivbcxit+1
          IF (vcxerr .GT. 0.05 .AND. ivbcxit .LT. 10)  go to 450
        END IF
!
        gi(j)    = erel0/ezero - (1.0_DP + 0.5_DP * taurat)*ge(j)
        IF (gi(j) .LT. zeroc)  gi(j) = zeroc
        gth(j)   = 1.0_DP - taurat*encap(j)
!
!       fionx allows testing of non-classical slowing down
!

        IF (fionx .EQ. zeroc)  go to 70
        ft    = fionx*gi(j)
        IF (ft .GT. ge(j))  ft = ge(j)
        gi(j) = gi(j) + ft
        ge(j) = ge(j) - ft
!
!       calculate bke and bki
!
   70   prat   = 1.0_DP

        IF (pzero(j) .NE. zeroc) prat = prel0/pzero(j)

        bke(j) = bkef_m(vrat,taurat,emzrat(j))*prat        ! this is just emzrat in onetwo

        bki(j) = prat - (1.0 + taurat)*bke(j)
  END DO      !  end loop over mesh points



  RETURN
!
  END       SUBROUTINE fast_ion_parms





      SUBROUTINE fast_ion_thermalization  (bke,encap,enfsav,ge,spbrfsav, & ! input
                        pprfsav,ppfsav,qfsav,sfsav,spfsav,taus,wfsav,    &
                        dtt,ibcur,nj,iangrot,                            &

                        enf,ppf,qf,sf,spf,taupf,                         & ! output         
                        tauppf,tauef,wf,pprf,spbrf)
! ----------------------------------------------------------------------
!     this is essentially slow2 in Onetwo 
!     This subroutine models fast ion slowing down.  The fast ions may
!     be beam ions from neutral beam heating or alpha particles from
!     fusion.  let the fast ion distribution function be given by
!         f(v,r,t) = g(v)*n(r,t).
!     This subroutine assumes that g(v) is fixed at all times at its
!     asymptotic value. Transients due to beam turn on and off are
!     accounted for only in n(r,t). That is, the fast ion density n(r,t)
!     at time t, is adjusted to produce the correct asymtotic source rate.
!     The change in g(v) over time to its asymtotic value is neglected.
!     for beam ions (not checked for fusion alphas) the quantities
!  INPUT
!     encap(1:nj)
!     bke(1:nj)
!     ge(1:nj)
!     taus(1:nj)
!     enfsav(1:nj)
!     sfsav(1:nj)
!     wfsav(1:nj)
!     qfsav(1:nj)
!     ppfsav(1:nj)
!     spfsav(1:nj)
!     pprfsav(1:nj)
!     spbrfsav(1:nj)
!     ibcur
!     iangrot
!     dtt
!     nj
!  OUTPUT :
!     taupf(1:nj)
!     tauppf(1:nj)
!     tauef(1:nj)
!     enf(1:nj)          fast ion density over rho grid (1/cm^3)
!     sf(1:nj)
!     qf(1:nj)
!     ppf(1:nj)
!     spf(1:nj)
!     spbrf(1:nj)
!     pprf(1:nj)
!
! ----------------------------------------------------------------------
!

     IMPLICIT  NONE
!
     REAL(DP) bke(*), encap(*), enfsav(*), ge(*), ppfsav(*), qfsav(*), &
                sfsav(*), spfsav(*), taus(*), wfsav(*)
     REAL(DP) enf(*), ppf(*), qf(*), sf(*), spf(*), taupf(*), &
                tauppf(*), tauef(*), wf(*), &
                pprf(*),spbrf(*),spbrfsav(*),pprfsav(*)
     REAL(DP)   dtt
     INTEGER(I4B) ibcur,iangrot,j,nj



!  begin loop over mesh points
!
  DO j=1,nj
!      calculate slowing down times
       taupf (j) = taus(j)*encap(j)
       tauppf(j) = taus(j) * ABS (bke(j))
       tauef (j) = 0.5_DP * taus(j)*ge(j)
!
!  calculate fast ion particle density and delayed particle source
!
      enf(j) = enfsav(j) * EXP (-dtt/taupf(j)) &
             + sfsav(j)*taupf(j)*(1.0_DP-EXP (-dtt/taupf(j)))
      sf (j) = enf(j)/taupf(j)
!
!  calculate fast ion energy density and delayed energy source
!
      IF(tauef(j) .LE. zeroc)THEN
         PRINT *,"tauef .le. izero at j =",j
         PRINT *,"taus(j) =",taus(j)
         PRINT *,"ge(j) =",ge(j)
      ENDIF
      wf(j) = wfsav(j) * EXP (-dtt/tauef(j)) &
            + qfsav(j)*tauef(j)*(1.0_DP-EXP (-dtt/tauef(j)))
      qf(j) = wf(j)/tauef(j)
!
!  calculate density and delayed source of fast ion parallel momentum
!
      IF (ibcur .NE. izero)THEN
         ppf(j) = ppfsav(j) * EXP (-dtt/tauppf(j)) &
             + spfsav(j)*tauppf(j)*(1.0_DP-EXP (-dtt/tauppf(j)))
         spf(j) = ppf(j)/tauppf(j)
      ENDIF

      IF (iangrot .NE. izero)THEN
!
! --- fast ion angular momentum density and delayed source
!
         pprf (j) = pprfsav(j) * EXP (-dtt/tauppf(j)) &
           + spbrfsav(j)*tauppf(j)*(1.0_DP-EXP (-dtt/tauppf(j)))
         spbrf(j) = pprf(j)/tauppf(j)
      ENDIF
!
    ENDDO   ! end loop over mesh points

      RETURN
!
  END SUBROUTINE fast_ion_thermalization


  SUBROUTINE fi_thermal_source(atw_beam,vionz,r2capi,rcap,     &
                              sb,qb,wb,pb0,fbe,fbi,fbth,spbr,bke,bki,ppb, &
                              ibslow,enbmin,        &
                              enb,kj,ke,kb,ic,jb,nj,external_beam_cur,codeid, &
                              r,eps,delta,kappa,rmajor,zeff,angrot,           &
                              sbeam,sscxl,qbbe,qbbi,qbeame_rot,          &
                              qbeami_rot,qbeame,qbeami,sprbeame,sprbeami,&
                              ssprcxl,enbeam,enbeams,wbeam,curbi,curbe,curbet)
!-----------------------------------------------------------------------------
! -- INPUT
! --    ibslow
! --    kj
! --    ic
! --    jb
! --    atw_beam
! --    vionz(1:nj)
! --    r2capi        <R**2>   ==> r2capi(1..nj)
! --    rcap          <R>      ==> rcap(1..nj)
! --    ke   declared second dim of pb0, etc, (=3, for full half,third energy)
! --    kb   declared third dim of pb0,etc. (number of beams)
! --    pb0(1:kj,1:ke,1:nbeams)
! --    sb(1:kj,1:ke,1:nbeams)
! --    spbr(1:kj,1:ke,1:nbeams)
! --    fbe(1:kj,1:ke,1:nbeams)
! --    fbi(1:kj,1:ke,1:nbeams)
! --    fbth(1:kj,1:ke,1:nbeams)
! --    qb(1:kj,1:ke,1:nbeams)
! --    bki(1:kj,ke,kb)
! --    bke(1:kj,ke,kb)
! --    enb(1:kj,ke,kb)
! --    wb(1:kj,ke,kb) 
! --    ppb(1:kj,ke,kb) 
! --    r(1:nj)
! --    zeff(1:nj)
! --    eps(1:nj)
! --    angrot
! --    external_beam_cur
! --    codeid
! --    kappa
! --    rmajor
! --    enbmin              minimum fast ion density !/cm^3
! -- OUTPUT
! --    sbeam(1:nj)    fast ions thermal paricle source #/(cm^3 sec)
! --    sscxl(1:nj)
! --    qbbe(1:nj,ke,kb)
! --    qbbi(1:nj,ke,kb)
! --    qbeame_rot(1:nj)
! --    qbeami_rot(1:nj)
! --    qbeame(1:nj)
! --    qbeami(1:nj)
! --    sprbeame(1:nj)
! --    sprbeami(1:nj)
! --    ssprcxl(1:nj)
! --    enbeam(1:nj)
! --    enbeams(1:nj)
! --    wbeam(1:nj)         fast ion stored energy density kev/cm3
! --    curbi(1:nj)
! --    curbet(1:nj)
! --    curbe(1:nj)
! --    delta
!-----------------------------------------------------------------------------

      IMPLICIT  NONE
      INTEGER(I4B) j,jb,ic,nj,kj,ke,kb,external_beam_cur,ibslow
      REAL(DP) erot,trot,atw_beam,prat,fscxl,fprscxl,delta,kappa,    &
               rmajor,enbmin
      REAL(DP) vionz(kj),r2capi(kj),rcap(kj),r(kj),eps(kj)
      REAL(DP) pb0(kj,ke,kb),sb(kj,ke,kb),qbbe(kj,ke,kb),qbbi(kj,ke,kb), &
               sbeam(kj),qbeame_rot(kj),qbeami_rot(kj),bke(kj,ke,kb),    &
               bki(kj,ke,kb),fbe(kj,ke,kb),fbi(kj,ke,kb),fbth(kj,ke,kb), &
               qbeame(kj),qbeami(kj),enbeam(kj),enbeams(kj),wbeam(kj),   &
               wb(kj,ke,kb),sscxl(kj),sprbeame(kj),sprbeami(kj),         &
               ssprcxl(kj),zeff(kj),enb(kj,ke,kb),curbi(kj),curbe(kj),   &
               curbet(kj),qb(kj,ke,kb),spbr(kj,ke,kb),ppb(1:kj,ke,kb),   &
               angrot(kj)

      CHARACTER(*) codeid


      DO  j=1,nj
         erot      = 0.5_DP * atw_beam*xmassp*vionz(j)**2
         trot      = atw_beam*xmassp*vionz(j)*r2capi(j)/rcap(j)
         prat      = 1.0_DP
         IF (pb0(j,ic,jb) .NE. zeroc) &
           prat    = 1.0_DP - atw_beam * xmassp * vionz(j) / pb0(j,ic,jb)

         ! --- fraction that doesn't thermalize is assumed lost due to
         ! ---  charge exchange . fbth is calculated in slow1 (as gth(j))
         ! --- and is taken as (1.0-(taus/taucx)*N ) where N is evaluated in
         ! --- function encapf and has the form
         ! ---     N = (taucx/taus)*(1-EXP (-tauf/taucx))
         ! --- which results when taucx is assumed independent of fast ion speed.
         ! --- Hence fscxl is simply
         !        (taus/taucx) * N = 1.0 - EXP (-tauf/taucx)             (HSJ)

         fscxl   = 1.0_DP  - fbth(j,ic,jb)
         fprscxl   = prat - bke (j,ic,jb) - bki(j,ic,jb)

         ! --- fraction (of fast ions) that do thermalize are a thermal ion source (HSJ):

         sbeam(j)  = sbeam(j) + (1.0_DP - fscxl) * sb(j,ic,jb)

         ! --- fast ion charge exchange rate(note that sscxl is not included
         ! --- as a particle source term. only the associated energy is
         ! --- accounted for (HSJ)

         sscxl(j)      = sscxl(j) + fscxl*sb(j,ic,jb)
         qbbe(j,ic,jb) = fbe(j,ic,jb)*qb(j,ic,jb)                  &
                    + bke(j,ic,jb)*spbr(j,ic,jb)*angrot(j)*1.0e-07 ! w/cm3
         qbbi(j,ic,jb) = fbi(j,ic,jb)*qb(j,ic,jb)                  &
                    + bki(j,ic,jb)*spbr(j,ic,jb)*angrot(j)*1.0e-07 &
                    + fbth(j,ic,jb)*sb(j,ic,jb)*erot*1.0e-07
         qbeame_rot(j) = qbeame_rot(j)                             &
                       + (bke(j,ic,jb)*spbr(j,ic,jb)*angrot(j)*1.0e-07) &
                       * 0.62415064e16
         qbeami_rot(j) = qbeami_rot(j)                                  &
                       + ( bki(j,ic,jb)*spbr(j,ic,jb)*angrot(j)*1.0e-07 &
                       + fbth(j,ic,jb)*sb(j,ic,jb)*erot*1.0e-07)        &
                       * 0.62415064e16

            qbeame(j)     = qbeame(j) + &
                              qbbe(j,ic,jb)*0.62415064e16 ! keV/(cm3 sec)
            qbeami(j)     = qbeami(j) + qbbi(j,ic,jb)*0.62415064e16

         sprbeame(j)   = sprbeame(j) + bke(j,ic,jb)*spbr(j,ic,jb)
         sprbeami(j)   = sprbeami(j) + bki(j,ic,jb)*spbr(j,ic,jb) &
                                     + fbth(j,ic,jb)*sb(j,ic,jb)*trot
         ssprcxl(j)    = ssprcxl(j)  + fprscxl*spbr(j,ic,jb) &
                                     + fscxl*sb(j,ic,jb)*trot


         IF (ibslow .NE. izero) THEN
            enb(j,ic,jb) = MAX(enb(j,ic,jb),enbmin)
            enbeam(j) = enbeam(j) + enb(j,ic,jb)
            enbeams(j)= enbeam(j)
            wbeam (j)  = wbeam(j) + kevpjou*wb(j,ic,jb)   !kev/cm3

 


            ! ----------------------------------------------------------------------
            !     BEAM-DRIVEN CURRENT CALCULATIONS
            ! ----------------------------------------------------------------------


            IF(external_beam_cur .EQ. 0)THEN
               curbi(j)  = curbi(j) + charge*ppb(j,ic,jb)   &
                                      /(2.99792458e9_DP*atw_beam*xmassp)
            ELSE
               curbi(j) = zeroc
            ENDIF

            !             Compute beam-driven electron current per rational fit to
            !             numerical results of D.F.H. Start and J.G. Cordey,
            !             Phys. Fluids, 23, 1477 (1980).

            curbe(j) = -curbi(j)/zeff(j)
            IF (codeid .EQ. 'onedee') THEN
               delta = r(j)/(SQRT (kappa)*rmajor)
            ELSE

               !                eps is advanced in time in rhomesh. HSJ

               delta = eps(j)
            END IF


            IF(external_beam_cur .EQ. 0)THEN
               curbet(j) = -curbe(j)*((1.55_DP+0.85_DP/zeff(j)) * SQRT (delta) &
                                     -(0.20_DP+1.55_DP/zeff(j)) *       delta)
            ELSE
               curbet(j) = 0.0 ! dont have this value from trasnp at present
            ENDIF           !so we just use total beam current
         ENDIF
      ENDDO   ! j loop over grid

      RETURN

   END SUBROUTINE fi_thermal_source



      REAL(DP) FUNCTION encapf_m (vcvo, tstcx)
! ----------------------------------------------------------------------
!  This routine evaluates the function N, which is related to the rate
!     at which fast ions slow down on electrons and thermal ions.  This
!     function is defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974).
!     It is assumed that taus/taucx is independent of the ion speed.
!     This allows calculation of Pcx,the probability against charge
!     exchange, as
!           Pcx = [(v0**3+vc**3)/(v**3+vc**3)]**(-taus/(3.0*taucx))
!     Subsequently the integral defining N can be evaluated with the
!     result that
!         N = (taucx/taus)*(1.0-EXP (-tauf/taus))
! --- input
!  tstcx            ratio of taus/taucx
!                   for fast alpha particle this ratio is taken as 0.0
!  vcvo             SQRT (Ec/E0)
! ------------------------------------------------------------------ HSJ
!
      IMPLICIT  NONE
      REAL(DP) v3,vcvo,tstcx,tfts

      v3   = vcvo**3
      tfts = LOG (1.0_DP + 1.0_DP/v3) / 3.0_DP
      IF (tstcx .LE. 0.01_DP) THEN
        encapf_m = tfts
      ELSE
        encapf_m = (1.0_DP - EXP (-tfts*tstcx)) / tstcx
      END IF

      RETURN

      END FUNCTION encapf_m



      REAL(DP)  FUNCTION bkef_m (vpar, tpar, emzpar)
! ----------------------------------------------------------------------
! This routine evaluates the function Ke, which is related to the
! transfer of parallel momentum from fast ions slowing down on electrons.
! The function is defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974).
! ----------------------------------------------------------------------

      IMPLICIT  NONE


      REAL(DP) vpar,tpar,emzpar,a,b,ep
      INTEGER(I4B) m4,m

!      REAL(DP),EXTERNAL  ::   bkefun_m  don't do this inside modules

          vcvo   = vpar
          tstcx  = tpar
          emzrat = emzpar

          a      = zeroc
          b      = 1.0_DP
          ep     = 0.1_DP
          m4     = 3

          bkef_m   = asimp_m (a, b, ep, m, m4,bkefun_m)


      RETURN

      END FUNCTION bkef_m


      REAL(DP) FUNCTION bkefun_m (y)
! ----------------------------------------------------------------------
! evaluate argument of integral defining Ke, the fraction of initial
! fast ion parallel momentum collisionally transferred to thermal electrons.
! see comments under function gefun_m HSJ
! ----------------------------------------------------------------------

      IMPLICIT NONE
      REAL(DP) y,v3,arg,alogarg,pcxlog,alogy3v3,alog3y,blog,    &
               bkeflog


      bkefun_m = zeroc
      IF (y .GT. zeroc) THEN
        v3       = vcvo**3
        arg      = (1.0 + v3)/(y**3+v3)
        alogarg  = LOG (arg)
        pcxlog   = -tstcx*alogarg*0.33333333333334_DP
        alogy3v3 = LOG (y**3+v3)
        alog3y   = 3.0_DP * LOG (y)
        blog     = (alog3y+alogarg)*emzrat*0.33333333333334_DP
        bkeflog  = alog3y+pcxlog+blog-alogy3v3
        IF (bkeflog .LT. -30.0_DP) THEN
          bkefun_m = zeroc
        ELSE
          bkefun_m = EXP (bkeflog)
        END IF
      END IF
      RETURN

      END FUNCTION bkefun_m



      REAL(DP)  FUNCTION gef_m (vpar, tpar)
! ----------------------------------------------------------------------
!  This routine evaluates the function Ge, which is related to the
!     transfer of energy from fast ions slowing down on electrons.  The
!     function is defined by Callen et al., IAEA Tokyo, Vol. I, 645
!     (1974).
!----------------------------------------------------------------------

      IMPLICIT  NONE
      INTEGER(I4B) m4,m
      REAL(DP) vpar,tpar,y,term2,term3,term4,rooty,a,b,ep
      REAL(DP) root13,fact4
!      REAL(DP), EXTERNAL  :: gefun_m don't do this inside modules



      DATA           root13 /0.577350_DP/, fact4 /0.302300_DP/



      IF (tpar .LE. 0.01_DP) THEN
        rooty = 1.0_DP / vpar
        y     = rooty**2
        term2 = LOG ((1.0_DP-rooty+y)/(1.0_DP + rooty)**2)/(6.0_DP * y)
        term3 = root13* ATAN (root13*(2.0_DP*rooty-1.0_DP))/y
        term4 = fact4/y
        gef_m   = 2.0_DP * (0.5_DP-term2-term3-term4)
      ELSE
        vcvo  = vpar
        tstcx = tpar
        a     = zeroc
        b     = 1.0_DP
        ep    = 0.1_DP
        m4    = 3
        gef_m   = asimp_m(a,b,ep,m,m4,gefun_m)
!        gef_m   = asimp_m(a,b,ep,m,m4,gefun_mm)
!        gef_m   = asimp_mmm(a,b,ep,m,m4,gefun_mmm)
      END IF
      RETURN

      END FUNCTION gef_m





      REAL(DP)  FUNCTION gefun_m (y)
! ----------------------------------------------------------------------
!  gefun_m is the argument of the Ge integral given by Callen.  the
!     substitution y = v/v0 has been made and the ratio of slowing down
!     time to charge exchange time, tstcx, is assumed to be independent
!     of speed.  the integral determing probability against charge
!     exchange can then be done analytically and leads to the expression
!     for pcx given here.  gefun_m is integrated using subroutine ASIMP.
!     for small densities (such as occur in large h mode gradients) tstcx
!     can become quite large (>1000) which causes numerical problems.
!     hence this routine was modified 4/26/89 by HSJ
! ----------------------------------------------------------------------


      IMPLICIT NONE
 
      REAL(DP) y,v3,arg,pcxlog,alogy3v3,gefunlog


      IF (y .LE. zeroc) THEN
        gefun_m = zeroc
      ELSE
        v3       = vcvo**3
        arg      = (1.0 + v3)/(y**3+v3)
        pcxlog   = -tstcx*LOG (arg)*0.33333333333334
        alogy3v3 = LOG (y**3+v3)
        gefunlog = Ln_2 + 4.0_DP * LOG (y) + pcxlog - alogy3v3
        IF (gefunlog .LT. -30.0_DP) THEN
          gefun_m = zeroc
        ELSE
          gefun_m = EXP (gefunlog)
        END IF
      END IF

      RETURN

      END FUNCTION gefun_m




      REAL(DP)  FUNCTION asimp_m (a1, b, ep, m, n, FUN)
!------------------------------------------------------------------------------
! --- author:  k. hillstrom (argonne national laboratory, chicago, illinois)
!------------------------------------------------------------------------------

      USE io_gcnmp,                                     ONLY : ncrt

      IMPLICIT NONE 



     REAL(DP),EXTERNAL :: FUN  ! no external statement needed  because FUN is dummy name
 
      REAL(DP)    a1, b, ep, a, eps, absar, est, fa, fm, fb, dx, sx
      REAL(DP)    f1,  est1, sum, daft, esum, tsum, da
      REAL(DP)    fmax, aest1, delta, aest,diff
      REAL(DP)    f2(30), fbp(30), est2(30), nrtr(30), aest2(30), ftst(3)


      INTEGER(I4B) m,n,i,lvl,l,n7

!     the parameter setup for the initial call

      n7 = ncrt
      IF (n .LE. 0) THEN
        WRITE (n7, '(a)')  ' **** ASIMP error return:  n .le. 0 ****'
        asimp_m = zeroc
        RETURN
      END IF

      IF (n .GT. 3) THEN
        WRITE (n7, '(a)')  ' **** ASIMP error return:  n .gt. 3 ****'
        asimp_m = zeroc
        RETURN
      END IF

      a       = a1
      eps     = ep*15.0
      esum    = zeroc
      tsum    = zeroc
      lvl     = 1
      da      = b-a
      fa      = FUN(a)
      fm      = FUN((a+b)*0.5)
      fb      = FUN(b)
      m       = 3
      fmax    = ABS (fa)
      ftst(1) = fmax
      ftst(2) = ABS (fm)
      ftst(3) = ABS (fb)
      DO 10 i=2,3
        IF (fmax .GE. ftst(i))  go to 10
        fmax = ftst(i)
   10 CONTINUE
      est   = (fa+4.0*fm+fb)*da/6.0
      absar = (ftst(1)+4.0*ftst(2)+ftst(3))*da/6.0
      aest  = absar

! 1 = recur

   20 dx         = da/(2.0**lvl)
      sx         = dx/6.0
      f1         = FUN(a+0.5*dx)
      f2(lvl)    = FUN(a+1.5*dx)
      est1       = sx*(fa+4.0*f1+fm)
      fbp(lvl)   = fb
      est2(lvl)  = sx*(fm+4.0*f2(lvl)+fb)
      sum        = est1+est2(lvl)
      ftst(1)    = ABS (f1)
      ftst(2)    = ABS (f2(lvl))
      ftst(3)    = ABS (fm)
      aest1      = sx*(ABS (fa)+4.0*ftst(1)+ftst(3))
      aest2(lvl) = sx*(ftst(3) +4.0*ftst(2) + ABS (fb))
      absar      = absar-aest+aest1+aest2(lvl)
      m          = m+2
      go to (60,30,70),n
   30 delta = absar
      go to 90
   60 delta = 1.0
      go to 90
   70 DO 80 i=1,2
        IF (fmax .GE. ftst(i))  go to 80
        fmax = ftst(i)
   80 CONTINUE
      delta = fmax
   90 diff  = ABS (est-sum)
      daft  = (est-sum)/15.0
      IF (diff-eps*delta) 110,110,100
  100 IF (lvl-30) 140,120,120
  110 IF (lvl- 1) 120,140,120

! 2 = up

  120 a    = a+2.0*dx
  130 lvl  = lvl-1
      esum = esum+daft
      l    = nrtr(lvl)
      tsum = tsum+sum
      go to (160, 170), l

! 11 = r1,12=r2

  140 nrtr(lvl) = 1
      est = est1
      aest = aest1
      fb = fm
      fm = f1
      eps = eps/2.0
  150 lvl = lvl+1
      go to 20
  160 nrtr(lvl) = 2
      fa = fb
      fm = f2(lvl)
      fb = fbp(lvl)
      est = est2(lvl)
      aest = aest2(lvl)
      go to 150
  170 eps = 2.0 * eps
      sum = zeroc
      IF (lvl-1) 180,180,130
  180 asimp_m = tsum
      a = ABS (esum)
      ep = diff/delta
      IF (a .GE. ep)  go to 190
      asimp_m = asimp_m - esum

  190 RETURN

      END FUNCTION asimp_m 




      REAL(DP) FUNCTION cxr_m (x)
! ----------------------------------------------------------------------
! this function calculates the charge exchange rate for hydrogen atoms
! interacting with protons in units of cm**3/s.
! x is in units of keV for 1.0e-3 .le. x .le. 100, the
! the formula is taken from the paper by r.l. freeman and e.m. jones
! clm-r 137 culham laboratory 1974.
! for x .lt. 1.0e-3, a rate coefficient derived from an analyti! average
! over the approximate cross section  sigma = 0.6937e-14*(1.0-0.155*LOG10
! (e/1ev))**2 is used.  this cross section is an approximation to that
! given by riviere, nuclear fusion 11,363(1971).
! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(DP) x,tc,dum

      IF (x .LT.   1.0e-3_DP)  go to 10
      IF (x .GT. 100._DP     )  go to 20
      tc  = LOG (x)+6.9077553
      dum = 0.2530205e-4-tc*0.8230751e-6
      tc  = -0.1841757e2+tc*(0.528295-tc*(0.2200477-tc*(0.9750192e-1-tc*   &
           (0.1749183e-1-tc*(0.4954298e-3+tc*(0.2174910e-3-tc*dum))))))
      tc  = EXP (tc)
      cxr_m = tc
      RETURN

   10 tc  = 0.50654_DP - 6.7316e-2_DP * LOG  (x)
      cxr_m =           3.3340e-7_DP * SQRT (x) * (tc*tc + 7.454e-3_DP)
      RETURN

   20 cxr_m = 0.0
      RETURN

      END  FUNCTION cxr_m




      REAL(DP) FUNCTION cxrv_m (x)
! ----------------------------------------------------------------------
!  NOTE THAT FUNCTION CXR USES THE MAXWELLIAN AVERAGE SIGMA V,
!       THIS FUNCTION ( ie CXRV ) USES SIGMA*V INSTEAD.
!       HSJ 6/5/96
! this FUNCTION calculates the charge exchange rate for hydrogen atoms
! interacting WITH protons in units of cm**3/s.
! x is in units of keV.
! for 1.0e-3 .LE. x .LE. 100, the
! the formula is taken from the paper by r.l. freeman and e.m. jones
! clm-r 137 culham laboratory 1974
! for x .LT. 1.0e-3, a rate coefficient derived from an analytic average
! over the approximate cross section  sigma = 0.6937e-14*(1.0-0.155*LOG10
! (e/1ev))**2 is used.  this cross section is an approximation to that
! given by riviere, nuclear fusion 11,363(1971).
! ----------------------------------------------------------------------

      IMPLICIT NONE
      REAL(DP) x,ev
      IF      (x .LT.   1.0e-3_DP) THEN
        PRINT *,'function CXRV_M: X is out of range: ',x
        lerrno = 355+ iomaxerr
        CALL terminate(lerrno,nlog)
      ELSE IF (x .GT. 100.0_DP  )  THEN
        cxrv_m = zeroc
      ELSE
        ev   = 1000.0_DP*x
        cxrv_m = 0.6937e-14_DP * (1.0_DP - 0.156_DP * LOG10 (ev))**2 /  &
             (1.0_DP + 0.1112e-14_DP * (ev**3.3))
      END IF
      RETURN

      END  FUNCTION cxrv_m


  END MODULE thermalization


