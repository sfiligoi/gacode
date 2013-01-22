      subroutine neut(taupm,srecom,sion,qione,qioni,qcx)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  Compute neutral source & reaction rates for charge exchange, 
c  electron ionization, and radiative recombination.
c  Note: ratef is a correction factor for elongated plasmas.
c
      implicit none
      include '../inc/data.m'
      include '../inc/ptor.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/input.m'
c
      logical reqsp
      integer icount, i, j, k, kion, nprim
      integer ineut(2), nneu, m1, in, iother, inrad, ibrem
      parameter(kion=5)  ! 3 primary, 2 impurity ions maximum
      real*8 zpi, volfac, sfarea, raneut, ratef, recrat, cxmix
      real*8 atw(kion), z(jmaxmm,kion),
     &       eionr(jmaxmm), eirate(jmaxmm),
     &       cexr(jmaxmm,2), cx12r(jmaxmm), recomr(jmaxmm)
      real*8 volsn(jmaxmm,2), srecom(jmaxmm,2)
      real*8 rm(jmaxmm), en(jmaxmm,kion), ene(jmaxmm), 
     &       enneu(jmaxmm,2), te(jmaxmm), ti(jmaxmm),
     &       tineu(jmaxmm)
      real*8 fluxn(2), taupm, twall, reflec, wrad, wion
      real*8 dn1(jmaxmm,2), dn2(jmaxmm,2), dnv(jmaxmm,2),
     &       wn1(jmaxmm,2), wn2(jmaxmm,2), wnv(jmaxmm,2),
     &       enn(jmaxmm,2), ennw(jmaxmm,2), tn(jmaxmm,2),
     &       tn_save(jmaxmm,2)
      real*8 sion(jmaxmm,2), sione(jmaxmm), qione(jmaxmm), 
     &       qioni(jmaxmm), scx(jmaxmm,2), qcx(jmaxmm),
     &       qrad(jmaxmm), sbcx(jmaxmm,2), qcxbl(jmaxmm)
c
      save icount, tn_save
c
      icount=icount+1
c
      raneut = dsqrt (kappa_d) * amin_d
      ratef = raneut / r_d(nj_d)
c
      in = 1
      nneu = 1
      reqsp = .true.
      inrad = 0
      ibrem = 1        ! don't include Bremsstrahlung radiation in qrad
      reflec = 1.D0    ! fraction of neutrals that reflect from wall
      twall = 1.D-2    ! temperature (keV) of recycled or reflected neutrals from wall
c     if (taupm.eq.0.) taupm = 0.35D0  ! ptcle confinement time
      zpi = atan2 ( 0.0D0, -1.0D0 )
      volfac = 4.D0*zpi**2.D0*rmajor_d*100.D0
      sfarea = volfac*hcap_d(jmaxm)*r_d(jmaxm)*100.D0
c
c... indicate whether species i has neutrals present
c
      ineut(1)=1
      ineut(2)=0
c
      nprim = 1    ! 1 primary ion for now
      do j=1,kion
        atw(j)=0.D0
      enddo
      atw(1)=2.D0  ! primary
      atw(2)=12.D0 ! impurity
c
c... set electron energy loss (keV) per ionization due to radiation
c    and ionization
c
      wrad = 0.006D0
      wion = 0.014D0
c
c... set charge parameters for primary ions
c    and assign densities (en(k,1)=thermal ion density
c    and en(k,2)=impurity density, ene(k)=electron density)
c    Note : only main ion density is transported
c    impurity and beam ion density taken from iterdb file
c
c     write(*,*) 'nj_d = ',nj_d
      write(*,*) 'jmaxm = ',jmaxm
      write(*,*) 'mxgrid = ',mxgrid
c     write(*,*) 'ngrid = ',ngrid
c      write(*,*) 'dinpro(ngrid+1) = ',dinpro(ngrid+1)
c      write(*,*) 'denpro(ngrid+1) = ',denpro(ngrid+1)
      write(*,*) 'ne_exp(mxgrid) = ',ne_exp(mxgrid)
c
      do k=1,jmaxm+1
         do i=1,nprim
           z(k,i) = 1.D0
           if (atw(i).ge.4.0) z(k,i) = 2.D0
           scx(k,i)=0.D0
           sbcx(k+1,i)=sbcx_exp(k,i)/1.D6
         enddo
         rm(k)=r_d(k)*100.D0
c         ene(k+1)=denpro(k)*1.D14                     ! n_e (1/cm^3)
c         en(k+1,1)=dinpro(k)*1.D14                    ! n_i
c         en(k+1,2)=en_d(k,2)*1.D-6                    ! n_Z
         ene(k+1)=ne_exp(k)*1.D14                      ! n_e (1/m^3)
         en(k+1,1)=ni_exp(k)*1.D14                     ! n_i
         en(k+1,2)=nz_exp(k)*1.D14                     ! n_Z
c         te(k+1)=Tepro(k)
c         ti(k+1)=Tipro(k)
         te(k+1)=te_exp(k)
         ti(k+1)=ti_exp(k)
         sione(k)=0.D0
         qione(k)=0.D0
         qioni(k)=0.D0
         qcx(k)=0.D0
         qrad(k)=0.D0
      enddo
      ene(1)=ene(2)
      en(1,1)=en(2,1)
      en(1,2)=en(2,2)
      te(1)=te(2)
      ti(1)=ti(2)
      sbcx(1,1)=sbcx(2,1)
      sbcx(1,2)=sbcx(2,2)
c
      do k=1,jmaxm+1
         write(*,22) k, r_d(k), ene(k), te(k), en(k,1), en(k,2)
      enddo
22    format('rec',2x,i2,2x,0p1f8.4,1p6e13.5)
c     write(*,*) 'rm(jmaxm) = ',rm(jmaxm)
c     write(*,*) 'hcap(jmaxm) = ',hcap_d(jmaxm)
c     write(*,*) 'sfarea = ',sfarea
c
c... get electron ionization and charge exchange <sigma-v>
c    reaction rate for neutrals (cm^3/s)
c
      call rate (atw,kion,jmaxmm,jmaxm,ratef,eionr,eirate,cexr,cx12r)
c
c... get radiative recombination <sigma-v> reaction rate of 
c    thermal ions and electrons into neutrals (cm^3/s)
c
      do k=1,jmaxm+1
        recomr(k) = ratef*recrat(te(k))*ene(k)
        write(*,22) k, r_d(k), ene(k), te(k), ratef, recomr(k)
      enddo
c
c... Compute neutral volume source, volsn, due to radiative recombination
c    and charge exchange of fast ions with species i.
c    Here, sbcx is source of thermal neutrals and source
c    of fast ions due to beam neutrals (#/(cm^3*s).
c    Compute qcx = ion power density due to neutral-ion charge exchange.
c
c    First piece computed below : ion charge exchange loss with beam neutrals
c    other piece not included : qcx(j) = qcx(j) - 1.5*tn(j,i)*fenn*sscxl(j)*ibcx
c    where fenn=1, ibcx=1, and sscxl is rate at which fast ions charge 
c    exchange with thermal neutrals to form thermal ion and fast neutral. 
c    This 2nd piece results in a GAIN for ion distribution due to fast ion 
c    thermal neutral charge exchange. To convert from keV/sec/cm*3 > watts/m**3, 
c    multiply qcx by 1.60217733e-10. Note: -qcx (W/m^3) in iterdb file.
c
      do k=1,jmaxm+1
        do i=1,nprim
          if (i.le.2) then
            srecom(k,i) = recomr(k)*z(k,i)*en(k,i)
            volsn (k,i) = -sbcx(k,i) + srecom(k,i)
            qcx(k)      = qcx(k) - 1.5*ti(k)*(sbcx(k,i)) 
            write(*,999) k,j,z(k,i),en(k,i),sbcx(k,i),
     &                   volsn(k,i),srecom(k,i),qcx(k),
     &                   -qcx_d(k)/1.60217733e-10
          endif
        enddo
      enddo
c
c... next, get fluxn
c    need gasflx(m,i) Flux of injected neutral gas for primary species i (1/cm^2-s)
c    nj -> jmaxm
c    r - > r_d*100
c    hcap - > hcap_d
c    ti -> Tipro
c
 2350 do 2351 i=1,2
 2351 if (ineut(i) .ne. 0) call copya (en(1,i),enneu(1,i),jmaxm)
      call copya (ti,tineu,jmaxm)
      do 2355 i=1,2
        if (ineut(i) .eq. 0)  go to 2355
        call trapv (r_d*100.D0,en(1,i),hcap_d,jmaxm+1,fluxn(i))
        do k=1,jmaxm+1
c         write(*,998) k, i, r_d(k)*100.D0, hcap_d(k), en(k,i), fluxn(i)
c         write(*,998) k, i, r_d(k)*100.D0, cexr(k,1), cexr(k,2), 
c    &                 cx12r(k), eionr(k)
        enddo
        fluxn(i) = fluxn(i)/(hcap_d(jmaxm+1)*r_d(jmaxm+1)*100.D0*taupm)
c       write(*,*) 'fluxn(i) = ',fluxn(i),hcap_d(jmaxm+1),
c    &             r_d(jmaxm+1)*100.D0, jmaxm
 2355 continue
c
 2360 m1 = 1
      if (nneu .eq. 1)  m1 = in
      call neucg2 (jmaxm+1,rm,ene,en(1,m1),en(1,2),ti,
     .          volsn(1,m1),volsn(1,2),eionr,cexr(1,m1),cexr(1,2),cx12r,
     .          nneu,reqsp,twall,atw(m1),atw(2),reflec,reflec,
     .          fluxn(m1),fluxn(2),
     .          dn1(1,m1),dn1(1,2),dnv(1,m1),dn2(1,1),dn2(1,2),dnv(1,2),
     .          wn1(1,m1),wn1(1,2),wnv(1,m1),wn2(1,1),wn2(1,2),wnv(1,2))
c
c     ineucg = ineucg + 1
c     ineu   = ineu+1
c
c     write(*,997) dn1(1,m1),dn1(1,2),dnv(1,m1),
c    .             dn2(1,1),dn2(1,2),dnv(1,2)
c     write(*,996) wn1(1,m1),wn1(1,2),wnv(1,m1),
c    .             wn2(1,1),wn2(1,2),wnv(1,2)
c
      call neuden(nneu,kion,jmaxmm,jmaxm+1,rm,hcap_d,volfac,volsn,
     &            sfarea, dn1,dn2,dnv,wn1,wn2,wnv,eirate,
     &            fluxn,cx12r,en,enn,ennw,tn)
c
c... save tn from 1st call which is used to compute qcxbl
c    from experimental qcx and then held fixed in time
      if (icount.eq.1) then
      do i=1,nprim
        do k=1,jmaxm+1
          tn_save(k,i)=tn(k,i)
        enddo
      enddo
      endif
c
      do i=1,nprim
        do k=1,jmaxm+1
c         write(*,995) k, i, r_d(k)*100.D0, en(k,i), enn(k,i),
c    &                 ennw(k,i), tn(k,i), cexr(k,1)
c    &         1.5*enn(k,i)*cexr(k,i)*en(k,i)*(Tipro(k)-tn(k,i))
c         write(*,999) k,j,rm(k),-1.5*Tipro(k)*(sbcx_exp(k,i)/1.D6),
c    &    -qcx_d(k+1)/1.60217733D-10
c    &    +1.5*Tipro(k)*(sbcx_exp(k,i)/1.D6)
c    &    -1.5*enn(k,i)*cexr(k,i)*en(k,i)*(Tipro(k)-tn(k,i)),
c    &         1.5*enn(k,i)*cexr(k,i)*en(k,i)*(Tipro(k)-tn(k,i)),
c    &         -qcx_d(k+1)/1.60217733D-10
        enddo
      enddo
c
 996  format('neut2-wn ',1p6e13.5)
 997  format('neut2-dn ',1p6e13.5)
 998  format('neut1 ',2i3,2x,0p1f8.4,1p6e13.5)
 995  format('neut2 ',2i3,2x,0p1f8.4,1p6e13.5)
 994  format('neut3 ',2i3,2x,0p1f8.4,1p6e12.4)
 989  format('qrad ',2i3,2x,0p1f8.4,1p6e12.4)
c
c... check for negative neutral densities
c
      do k=1,jmaxm+1
         if (nneu .eq. 1) then
            if (enn(k,in) .lt. 0.0)
     .      write (6, *) " ERROR: enn(k,in),k =", enn(k,in), k
         else
            if (enn(k,1) .lt. 0.0 .or. enn(k,2) .lt. 0.0)
     .      write (6, *) " ERROR: enn(k,1), enn(k,2), k =",
     .                            enn(k,1), enn(k,2), k
         end if
      end do
c
c... compute particle and energy sources due to recombination of ions
c    and ionization and charge exchange of neutrals
c    First, include beam electron source #/(m^3*s) in electron ptcle source
c    Here, qrad due to recombination and ionization should be added to
c    Bremsstrahlung and cyclotron losses.
c
      do k=1,jmaxm+1
        sione(k)=sbion_d(k)*1.D-6
      enddo
c
      do 2430 k=1,jmaxm+1
      do 2430 i=1,nprim
      if (i .gt. 2)  go to 2430
      sione(k) = sione(k) - z(k,i)*srecom(k,i)
      qione(k) = qione(k) + 1.5D0*z(k,i)*srecom(k,i)*te(k)
      qioni(k) = qioni(k) - 1.5D0*srecom(k,i)*ti(k)
      if (inrad .eq. 0) qrad(k) = qrad(k) + z(k,i)*srecom(k,i)*wion
      if (ibrem .ne. 0) qrad(k) = qrad(k) + 0.62415064D16*ene(k)*
     &          ( 5.34D-31*dsqrt(te(k))*en(k,i)*z(k,i)**2 )
      srecom(k,i) = -srecom(k,i)
      if (ineut(i) .eq. 0)  go to 2430
      sion(k,i) = enn(k,i)*eirate(k)
      sione(k) = sione(k) + z(k,i)*sion(k,i)
      qione(k) = qione(k) + z(k,i)*sion(k,i)*wion
      qioni(k) = qioni(k) + 1.5D0*sion(k,i)*tn(k,i)
c     write(*,994) k, i, rm(k), sion(k,i), sione(k),
c    &             srecom(k,i), qione(k), qioni(k)
      qcx(k) = qcx(k) + 1.5D0*enn(k,i)*cexr(k,i)*
     &         en(k,i)*(ti(k)-tn(k,i))
      qcxbl(k)=-qcx_d(k)/1.60217733D-10+
     &         1.5D0*ti_d(k)*sbcx(k,i)-
     &         1.5D0*en_d(k,i)*cexr(k,i)*enn_d(k,i)*
     &         (ti_d(k)-tn(k,i))/1.D12
      qcx(k)=qcx(k)+qcxbl(k)
c     write(*,994) k, i, tn(k,i), cexr(k,i), qcx(k)
c     write(*,994) k, i, rm(k), tn_save(k,i),
c    &  -1.5D0*ti(k)*sbcx(k,i), qcxbl(k),
c    &  1.5D0*enn(k,i)*cexr(k,i)*en(k,i)*(ti(k)-tn(k,i)),
c    &  -qcx_d(k)/1.60217733D-10, qcx(k)
      if (inrad .eq. 0)  qrad(k) = qrad(k) + z(k,i)*sion(k,i)*wrad
c     write(*,989) k, i, r_d(k)*100.D0, z(k,i)*srecom(k,i)*wion,
c    &      z(k,i)*sion(k,i)*wrad,
c    &      0.62415064D16*ene(k)*5.34D-31*dsqrt(te(k))*
c    &      en(k,i)*z(k,i)**2,qrad(k)
      if (nneu .eq. 1)  go to 2430
      iother = 3-i
      cxmix = en(k,i)*enn(k,iother)*cx12r(k)
c     write(*,*) 'cxmix = ',cxmix
      scx(k,i) = scx(k,i) - cxmix
      scx(k,iother) = scx(k,iother) + cxmix
      qcx(k) = qcx(k) + 1.5*cxmix*(ti(k)-tn(k,iother))
 2430 continue
c
c... for volsn:
c
c    compute neutral volume source due to radiative recombination
c    Here, need volsn which is the volume source of neutrals,
c    recomr, charge for primary ions z(j,i), and the primary ion 
c    densities en(j,i) where j=1,2..nj,i=1,2..nprim (#/cm**3)
c
c    In cray309, we have
c
c     fast beam ion sources due to charge exchange with species i, sbcx
c     call bproc (hicm, kb, ke, kj, nbeams, nj, sbsav, sbcx, sbion)
c
c     do i=1,2
c       do j=1,nj
c         volsn(j,i) = sbcx(j,i)            ! thermal neutral and fast..
c       end do
c     end do
c
c     Here, sbcx is the sink due to cx with beam neutrals (#/(m^3*s), species: d
c     and is written to iterdb file from cray204 as
c629                if (iterdsc .ne. 0)
c630      .         write  (niterdb, 4158) namep(jj)
c631  4158          format ('*  sbcx : sink due to cx with beam neut.,',
c632      .                 ' #/(meter**3*second), species: ', a)
c633                write  (niterdb, 10) (sbcx(j,jj)*1.0e+6, j=1,nj)
c
c     This is currently read from iterdb file from readxp.f and
c     is saved in the variable sbcx_d.
c
c
c
c     do j=1,njs
c       recomr = ratef*recrat(te(j))*ene(j)
c       do i=1,nprim
c         if (i .le. 2) then
c           srecom(j,i) = recomr*z(j,i)*en(j,i)
c           volsn (j,i) = volsn(j,i) + srecom(j,i)
c         end if
c       end do
c     end do
c
c... diagnostic printout
c
c     do k=1,ngrid
c       write(*,100) k, Tepro(k), Tipro(k), denpro(k)*1.D14,
c    &               eira(Tepro(k),atw(1)), cxra(Tipro(k),atw(1))
c       write(*,100) k, denpro(k)*1.D14, Tepro(k), ratef, recomr(k)
c     enddo
c
c... from cray309.f in ONETWO :
c
c    enn(j,i) = neutral density for species i (#/cm^3)
c    tn(j,i) = temperature of neutral density for species i
c    qione = electron power density due to recombination and 
c            impact ionization (keV/sec/cm*3)
c    qioni = ion power density due to recombination and
c            impact ionization (keV/sec/cm*3)
c    qcx = ion power density due to neutral-ion charge 
c            exchange (keV/sec/cm*3)
c    Here, multiply qione, etc by 1.60217733e-10 to go 
c    from (keV/sec/cm*3) to (W/m**3). The sign conventions in the
c    iterdb files are :
c       qione - negative
c       qioni - positive
c       qcx - negative
c
c     do j=1,ngrid
c       sion(j,i) = enn(j,i)*eirate(j)
c       sione(j) = sione(j) + z(j,i)*sion(j,i)
c       qione(j) = qione(j) + z(j,i)*sion(j,i)*wion
c       qioni(j) = qioni(j) + sion(j,i)*1.5*tn(j,i)
c       qcx(j) = qcx(j) + 1.5*enn(j,i)*cexr(j,i)*en(j,i)*(ti(j)-tn(j,i))
c     enddo
c
100   format(2x,i2,2x,1p6e12.4)
999   format('rec2 ',2i3,2x,0p1f8.4,1p6e13.5)
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine rate (atw,kion,jm,jmaxm,ratef,eionr,eirate,cexr,cx12r)
c
c     RATE computes the neutral reaction rate coefficients <sigma-v>
c     for electron ionization and charge exchange. The units for
c     the routines eira and cxra are CGS and are in cm^3/s. 
c     So, 1.e6 conversion factor added to put rates in MKS units.
c     Within XPTOR, denpro is the electron density=0.10*ne_exp
c     where ne_exp is in units of 10^19 m^-3. 
c     Note: ratef is the correction factor for elongated plasmas.
c
c     Only one neutral species for now ...
c
      implicit none
      include '../inc/data.m'
      include '../inc/ptor.m'
c
      integer k, kion, jm, jmaxm
c     parameter(kion=5)  ! 3 primary, 2 impurity ions
      real*8 ratef, eira, cxra
      real*8 atw(kion), eionr(jm), eirate(jm),
     &       cexr(jm,2), cx12r(jm),
     &       te(jm), ti(jm), ene(jm)
c
      do k=1,jmaxm+1
        eionr(k)  = 0.D0
        eirate(k) = 0.D0
        ene(k+1)=denpro(k)*1.D14  ! n_e
        te(k+1)=Tepro(k)          ! Te
        ti(k+1)=Tipro(k)          ! Ti
      enddo
      ene(1)=ene(2)
      te(1)=te(2)
      ti(1)=ti(2)
c
      do k=1,jmaxm+1
        eionr(k)    = ratef*eira(te(k),atw(1))
        eirate(k)   = eionr(k)*ene(k)
        cexr(k,1)   = ratef*cxra(ti(k),atw(1))
        cexr(k,2)   = 0.D0
        cx12r(k)    = 0.D0
      enddo
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function eira (te, atw)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c     EIRA computes the electron ionization rate coefficient
c     <sigma-v> for any primary ion: hydrogenic and helium ions
c     have to be handled separately.
c ----------------------------------------------------------------------
c
      if (atw .gt. 3.D0)  go to 10
c
c     hydrogenic
c
      eira = eir(te)
      return
c
c     helium
c
   10 eira = eirhe(te)
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function eir (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c This function calculates the ionization rate of atomic hydrogen by ele
c impact in units of cm**3/s.  x is in units of keV.
c The formula is taken from the paper by R.L. Freeman and E.M. Jones
c clm-r 137 culham laboratory, 1974.
c ----------------------------------------------------------------------
c
      if (x .lt. .001D0 .or. x .gt. 100.D0)  goto 20
      ta = dlog (x)+6.9077553D0
      ta = -0.3173850D2+ta*(0.1143818D2-ta*
     &     (0.3833998D1-ta*(0.7046692D0-ta
     &     *(0.7431486D-1-ta*(0.4153749D-2-ta*0.9486967D-4)))))
      ta = dexp (ta)
      eir = ta
      return
   20 eir = 0.D0
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function eirhe (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c This function calculates the ionization rate of atomic helium by elect
c impact in units of cm**3/s.  x is in units of keV
c The formula is taken from the paper by R.L. Freeman and E.M. Jones
c clm-r 137 Culham Laboratory, 1974.
c ----------------------------------------------------------------------
c
      if (x .lt.  0.025)  xx = dlog(0.025D0) + 6.9077553D0
      if (x .gt. 10.0  )  xx = dlog(10.D0) + 6.9077553D0
      if (x .ge. .025D0 .and. x .le. 10.D0) xx = dlog(x) + 6.9077553D0
      xx = -0.4450917D2 + xx*(0.2442988D2 + xx*(-0.1025714D2
     &     + xx*(0.2470931D1 + xx*(-0.3426362D0 + xx*(0.2505100D-1
     &     - xx*0.7438675D-3)))))
      eirhe = dexp(xx)
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function cxra (ti, atw)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c     CXRA returns the charge exchange rate coefficient <sigma-v>
c     for any primary ion: hydrogenic and helium ions have to be
c     handled separately.
c ----------------------------------------------------------------------
c
      if (atw .gt. 3.D0)  go to 10
c
c     hydrogenic
c
      cxra = cxr(ti/atw)
      return
c
c     helium
c
   10 cxra = cxrhe(ti)
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function cxr (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c ----------------------------------------------------------------------
c This function calculates the charge exchange rate for hydrogen atoms
c interacting with protons in units of cm**3/s.
c x is in units of keV for 1.0e-3 .le. x .le. 100,
c the formula is taken from the paper by R.L. Freeman and E.M. Jones
c clm-r 137 culham laboratory 1974.
c For x .lt. 1.0e-3, a rate coefficient derived from an analytic average
c over the approximate cross section  sigma = 0.6937e-14*(1.0-0.155*LOG10
c (e/1ev))**2 is used.  This cross section is an approximation to that
c given by riviere, Nuclear Fusion 11,363(1971).
c ----------------------------------------------------------------------
c
      if (x .lt. 1.D-3)  goto 10
      if (x .gt. 100.D0) goto 20
      tc  = LOG (x)+6.9077553D0
      dum = 0.2530205D-4-tc*0.8230751D-6
      tc  = -0.1841757D2+tc*(0.528295D0-tc*
     &      (0.2200477D0-tc*(0.9750192D-1-tc*
     &      (0.1749183D-1-tc*
     &      (0.4954298D-3+tc*(0.2174910D-3-tc*dum))))))
      tc  = dexp(tc)
      cxr = tc
      return
c
   10 tc  = 0.50654D0 - 6.7316D-2 * dlog(x)
      cxr = 3.3340D-7 * dsqrt(x) * (tc*tc + 7.454D-3)
      return
c
   20 cxr = 0.D0
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function cxrhe (tkev)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c This function calculates the double charge exchange
c rate for helium as a function of temperature in keV.
c The equations were derived by performing
c a least squares curve fit on data provided by:
c     Oak Ridge National Labratory
c     atomic data for controlled fusion, ornl-dwg 76-6797
c
c vel(thermal)  = sqr(2t/mass(he))
c 5.4676655e-11 = sqr(2 /mass(he))
c t = temperature in ergs
c mass is in grams
c ----------------------------------------------------------------------
c
      vth = 5.4676655D-11 * dsqrt(tkev/1.6D-12)
      atkev = dlog10(tkev)
      if (tkev .gt. 1.587D3) ratehe = 0.D0
      if (tkev .lt. 0.05D0)  ratehe = 4.D-16
      if (tkev .le. 90.D0 .and. tkev .ge. 0.05D0)
     &  ratehe = vth*(10.D0**(-15.6868D0+atkev*
     &           (-0.192101D0+atkev*2.21944D-2)))
      if (tkev .gt. 90.D0 .and. tkev .le. 1.587D3)
     &  ratehe = vth*(10.D0**(-24.9253D0+atkev*
     &           (9.07495D0-atkev*2.31366D0)))
      cxrhe  = ratehe
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function recrat (te)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     RECRAT computes the local rate of radiative recombination
c     of thermal ions and electrons into neutrals (#/cm**3-s).
c     The <sigma-v> reaction rate coefficient has been fit
c     by a rational polynomial to the formula in Gordeev et. al.
c     [JETP Lett.,V.25,No.4,p.204,1977] to better than 0.1%.
c     Multiply this rate by Z of the ions in the calling routine.
c     The limits on Te are as specified in Gordeev et.al.
c     Te is in keV.
c
      t = te
      if (t .lt. 1.D-3) t = 1.D-3
      if (t .gt. 25.D0) t = 25.D0
      t = dlog10(t)
      recrat = (-15.42043D0+t*(-6.215405D0+t*
     &         (-1.559212D0-t*0.8086629D-1)))/
     &         (1.D0 + t*(0.3166686D0+t*(0.6894623D-1-t*0.1031934D-2)))
      recrat = 10.D0**recrat
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine copya (old_array, new_array, number)
c
      implicit none
c
c --- subroutine COPYA copies the array OLD_ARRAY into the array NEW_ARRAY
c
      integer  index, number
      real*8   old_array(*), new_array(*)
c
      do index=1,number
        new_array(index) = old_array(index)
      end do
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine trapv (r, y, fact, nj, xint)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     this subroutine integrates fact(r)*r*y(r) with respect to r from
c     zero to rminor.  the trapezoidal rule is used.
c
      dimension  r(*), fact(*), y(*)
c
      xint = 0.0
      do 10 j=2,nj
   10 xint = xint+0.5*(r(j-1)*fact(j-1)*y(j-1)+r(j)*y(j)
     &      *fact(j))*(r(j)-r(j-1))
      return
c
      end
