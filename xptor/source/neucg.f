      subroutine neucg2 (nx, rext, dene, den1, den2, ti, vols1, vols2,
     .                   eionr, cx1r, cx2r, cx12r, nspec, rxeq, twall,
     .                   atw1, atw2, rd1, rd2, flux1, flux2, dn11, dn12,
     .                   dn1v, dn21, dn22, dn2v, wn11, wn12, wn1v, wn21,
     .                   wn22, wn2v)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     NEUCG2
c
c     Two Species Neutral Transport Code
c
c     Given the plasma ion density and temperature profiles, and
c     the recycling plasma flux, NEUCG2 calculates the steady-state
c     density and temperature (i.e. 2/3 average energy) profiles for
c     one or two neutral species.  The treatment includes electron
c     ionization and charge exchange reactions, but not ion impact
c     ionization.  The boundary condition at the plasma surface is
c     a diffuse reflection at the wall temperature.  We further
c     assume that the charge exchange rate coefficient <sigma-v>
c     is independent of the relative particle velocity.
c
c        Inputs
c
c        nx     - Number of points in external radial grid
c        rext   - External radial grid points (cm)
c        dene   - Electron density (1/cm**3)
c        den1   - Ion species 1 density (1/cm**3)
c        den2   - Ion species 2 density (1/cm**3)
c        ti     - Ion temperature (both species) (keV)
c        vols1  - Neutral volume source for species 1 (1/cm**3-sec)
c                 (recombination, beam deposition)
c        vols2  - Same for species 2
c        eionr  - Electron impact ionization rate coefficient
c                 <sigma-v> (cm**3/sec)
c        cx1r   - Charge exchange rate coefficient for species 1
c                 <sigma-v> (cm**3/sec)  (same species reactions only)
c        cx2r   - Same for species 2
c        cx12r  - Charge exchange rate coefficient for mixed species
c                 reactions (both ways) <sigma-v> (cm**3/sec)
c        nspec  - Number of neutral species (1 or 2)
c                 In 1 species mode, only species 1 quantities are used,
c                 and mixed species charge exchange reactions are ignored.
c        rxeq   - Logical flag: true= rext is equally spaced,
c                 false= rext is unevenly spaced.  If false, an equally
c                 spaced grid with nx points is generated internally to
c                 reduce execution time. The calling procedure is the
c                 same in all other respects.
c        twall  - Temperature of recycled plasma and reflected neutrals
c                 from the wall (keV)
c        atw1   - Atomic weight of species 1 (h = 1, d=2, t=3)
c        atw2   - Same for species 2
c        rd1    - Fraction of neutrals of species 1 that reflect back
c                 from the wall (0-1, where 1 conserves particles)
c        rd2    - Same for species 2
c        flux1  - Flux of recycled plasma of species 1 at surface
c                 (1/cm**2-sec)
c        flux2  - Same for species 2
c
c                            Outputs
c
c        dn11   - Neutral density of species 1 due to flux1 (1/cm**3)
c                 (proportional to flux1)
c        dn12   - Neutral density of species 1 due to flux2 (1/cm**3)
c                 (proportional to flux2 - due to mixed charge exchanges)
c        dn1v   - Neutral density of species 1 due to volume source
c                 (1/cm**3)
c        dn21   - Neutral density of species 2 due to flux1
c        dn22   - Neutral density of species 2 due to flux2
c        dn2v   - Neutral density of species 2 due to volume source
c        wn11   - Temperature (2/3 average energy) of neutrals dn11 (keV)
c        wn12,wn1v,wn21,wn22,wn2v - Similarly
c
c        Note that the total neutral density of species 1 = dn11+dn12+dn1v
c        and the total neutral density of species 2 = dn21+dn22+dn2v.
c        All outputs are computed at external grid points (rext).
c
c                           Parameters
c
c        nr     - Number of points in radial quadrature (fixed by data)
c        nt     - Number of points in theta quadrature (fixed by data)
c        nxmax  - Maximum allowed value of nx (recompile to change).
c                 NOTE THAT nxmax APPEARS IN MANY MODULES.
c
c                             Notes
c
c          Throughout the program, variable names with w are usually
c        energy related, p often represents rho, 1 or 2 names species,
c        o- (over) denotes reciprocal (saves on division), variables
c        in /nudiag/ ending with -s are stored for use by Nuspec, and
c        temp is for temporary.
c          Where appropriate, do loop index i represents the fixed
c        radius r=r(i), do loop index j represents the running
c        integration variable rho=p=r(j), and do loop index k
c        represents the base point in the theta integral.
c          In two-species mode, operations are done in parallel, where
c        possible, to save execution time.
c
c        NEUCG2 is a re-write of Keith Burrell's NEUCG, as modified
c        by H. Howe.  Neucg is described in Burrell, K., Journal of
c        Comp. Physics, vol. 27, no. 1, p. 88, April 1978.
c                                    - Bob Stockdale,  April 1982.
c ----------------------------------------------------------------------
c
      include '../inc/neucg.m'
c
      dimension rext(*),dene(*),den1(*),den2(*),ti(*)
      dimension vols1(*),vols2(*),eionr(*),cx1r(*),cx2r(*),cx12r(*)
      logical   rxeq, rxeqs
      dimension dn11(*),dn12(*),dn1v(*),dn21(*),dn22(*),dn2v(*)
      dimension wn11(*),wn12(*),wn1v(*),wn21(*),wn22(*),wn2v(*)
c
      dimension tii(nr)
      dimension pn11(nr),pn12(nr),pn1v(nr),pn21(nr),pn22(nr),pn2v(nr)
      dimension pw11(nr),pw12(nr),pw1v(nr),pw21(nr),pw22(nr),pw2v(nr)
      dimension fn11(nr),fn12(nr),fn1v(nr),fn21(nr),fn22(nr),fn2v(nr)
      dimension fw11(nr),fw12(nr),fw1v(nr),fw21(nr),fw22(nr),fw2v(nr)
      dimension temp(nxmax),temp1(nxmax)
c
      logical         twospc, onespc
      common /nuspcs/ twospc, onespc
      common /nuquad/ r(nr),wr(nr),rwr(nr),sint(nt),wt(nt)
      common /nuatr0/ atrr01(nr,nt),atar01(nr,nt),
     .                atrr02(nr,nt),atar02(nr,nt)
      common /nuatbl/ odelr, ra(nxmax), a1(nxmax), a2(nxmax), nra, nram1
      common /nuovth/ ovth1(nr),ovth2(nr)
      common /nukern/ fb12(nr),fb21(nr),k11(nr,nr),k11w(nr,nr),
     .                                  k22(nr,nr),k22w(nr,nr)
      real*8          k11,k11w,k22,k22w
      common /nukint/ ki1(nr),ki1w(nr),ki2(nr),ki2w(nr)
      real*8          ki1,ki1w,ki2,ki2w
      common /nuhphi/ h1(nr),h1w(nr),phi1(nr),h2(nr),h2w(nr),phi2(nr)
      common /nubvth/ b1v(nr),b2v(nr),vth1(nr),vth2(nr)
      common /nuflxi/ flxh1,flxh2,flxi1(nr),flxi2(nr)
      common /nudiag/ alph1,alph2,alphv,beta1,beta2,betav,
     .                tws,g1s,g2s,atw1s,atw2s,rxs(nxmax),
     .                tis(nxmax),b11s(nxmax),fb12s(nxmax),vn1s(nxmax),
     .                b22s(nxmax),fb21s(nxmax),vn2s(nxmax),
     .                rxeqs
      common /nuvolsrc/vn1(nr),vn2(nr)
      data pi/3.141592654/,sqrpi/1.772453851/,cmskev/4.376978e7/
c
c                           Initialization
c
c          Set up radial and theta quadrature grids and weights
c
      if (   nx .gt. nxmax)  call nuerr(1)
      if (flux1 .eq. 0.0  )  call nuerr(5)
      rminor = rext(nx)
      rm2 = rminor**2
      call nugrid(rminor)
      twospc = .true.
      if (nspec .ne. 2)  twospc = .false.
      onespc = .not.twospc
c     do i=1,nx
c       write(*,111) i, rext(i), dene(i), den1(i), den2(i), ti(i)
c       write(*,111) i, rext(i), eionr(i), cx1r(i), cx2r(i), cx12r(i)
c     enddo
 111  format('neucg ',i3,2x,0p1f8.4,1p6e13.5)
c
c          Set up table of A1 and A2: the neutral loss rates
c
      nra   = nx
      nram1 = nra - 1
c
      do 110 i=1,nra
        ra(i) = rext(i)
        if (onespc)  go to 105
        a1(i) = dene(i)*eionr(i) + den1(i)*cx1r(i) + den2(i)*cx12r(i)
        a2(i) = dene(i)*eionr(i) + den2(i)*cx2r(i) + den1(i)*cx12r(i)
        go to 110
  105   a1(i) = dene(i)*eionr(i) + den1(i)*cx1r(i)
  110 continue
c
c     convert to an equally-spaced grid, if not already
c
      if (.not. rxeq)  call aspace(temp,temp1)
      rxeqs = rxeq
      odelr = 1.0 / ra(2)
c
c     interpolate quantities to internal grid
c
      call nuint (nx,rext,ti,nr,r,tii)
      call nuint (nx,rext,vols1,nr,r,vn1)
      if (twospc)  call nuint (nx,rext,vols2,nr,r,vn2)
      do 115 i=1,nx
 115  b11s(i) = den1(i)*cx1r(i)
      call nuint (nx,rext,b11s,nr,r,b1v)
c
      if (onespc)  go to 135
      do 120 i=1,nx
 120  b22s(i) = den2(i)*cx2r(i)
      call nuint (nx,rext,b22s,nr,r,b2v)
c
      do i=1,nx
        fb12s(i) = cx12r(i)/cx1r(i)
      end do
      call nuint (nx,rext,fb12s,nr,r,fb12)
c
      do i=1,nx
        fb21s(i) = cx12r(i)/cx2r(i)
      end do
      call nuint (nx, rext, fb21s, nr, r, fb21)
c
  135 do i=1,nr
        vth1(i) = cmskev * SQRT (tii(i)/atw1)
        ovth1(i) = 1.0 / vth1(i)
        vn1(i) = vn1(i)/b1v(i)
        b1v(i) = b1v(i)/(sqrpi*vth1(i))
        if (.not. onespc) then
          vth2(i) = cmskev * SQRT (tii(i)/atw2)
          ovth2(i) = 1.0 / vth2(i)
          vn2(i) = vn2(i)/b2v(i)
          b2v(i) = b2v(i)/(sqrpi*vth2(i))
        end if
      end do
c
      do i=1,nx
        rxs(i) = rext(i)
        tis(i) = ti(i)
        vn1s(i) = vols1(i)/b11s(i)
        if (twospc) vn2s(i) = vols2(i)/b22s(i)
      end do
c
      vthw1 = cmskev * SQRT (twall/atw1)
      swall1 = 2.0 * flux1/vthw1
      tws = twall
      g1s = swall1/(pi*vthw1**3)
      atw1s = atw1
      if (onespc)  go to 150
c
      vthw2 = cmskev * SQRT (twall/atw2)
      flux2i = flux2
      s2 = 1.0
      if (flux2 .eq. 0.0) s2 = 0.0
      if (flux2 .eq. 0.0) flux2i = 1.0
      swall2 = 2.0 * flux2i/vthw2
      g2s = swall2/(pi*vthw2**3)
      atw2s = atw2
c
c          Initialization completed.  We now have:
c
c        tii    - Ion temperature (keV)
c        vth1   - Thermal velocity for species 1 (cm/sec)
c        vth2   - Same for species 2
c        b1v    - Factor for species 1 kernel: B11/(sqrpi*vth1)
c                 (B11 is the species 1 birth rate due to same species
c                 charge exchange)
c        b2v    - Same for species 2
c        fb12   - Mixed species charge exchange conversion factor for
c                 species 1 kernel: B12/B11 (converts kernel 1 to
c                 compute birth rate of species 1 neutrals due to
c                 mixed species reactions)
c        fb21   - Similarly for species 2
c        vn1    - Neutral volume source conversion factor for species 1:
c                 vols1/B11
c                 (converts kernel to compute neutral volume source)
c        vn2    - Same for species 2
c        vthw1  - Thermal velocity of species 1 neutrals from wall source
c                 (recycled plasma and reflected neutrals)  (cm/sec)
c        vthw2  - Same for species 2
c        swall1 - Normalization of Maxwellian wall source of species 1 that
c                 produces a flux = flux1
c        swall2 - Same for species 2
c        s2     - Flag indicating whether flux2 = 0: When flux2=0 (s2=0)
c                 we use an arbitrary value for flux2 in order
c                 to compute the profile of reflected species 2 neutrals.
c                 Later we discard the terms for species 2 originating
c                 from the recycled plasma flux.
c        flux2i - Internal value of flux2 used (as described above)
c
c          These quantities have been saved for the neutral spectrum
c          routine (Nuspec) in /nudiag/:
c          rxeqs,tws,g1s,atw1s,g2s,atw2s,rxs,tis,b11s,b22s,vn1s,vn2s,
c          fb12s,fb21s.  Note that arrays are at external grid points.
c          Also passed to Nuspec is the A function table in /nuatbl/.
c          Later the reflection factors (alph1, etc.) will also be
c          saved in /nudiag/.
c
c ----------------------------------------------------------------------
c
c          Set up kernels and singular part of integrals
c
  150 call kernel
      call kintegral
c
c        The theta integral of each K(p,r) passes through a singulariy
c        at p = r.  This is handled by adding and subtracting the
c        transposed kernel K(r,p) (see write-up).  The total contribution
c        to the integral at the base point p = r is computed here.  This
c        is stored, conveniently, in the now empty diagonal of the K's.
c        (Note that at r = 0 the integral is not singular, and so the
c        quadrature is straightforward.  The special (and simple) term
c        at p = r=0 was computed by subroutine KINTEGRAL.)
c
c          Sum up the quadrature on the transposed kernel, and add
c          the singular part.  We use K(r = r(i),p=r(j)) = K(i,j).
c
      k11 (1,1) = ki1 (1)
      k11w(1,1) = ki1w(1)
      if (twospc) k22(1,1) = ki2(1)
      if (twospc) k22w(1,1) = ki2w(1)
c
      do i=2,nr
        sum = 0.0
        do j=1,nr
          sum = sum + rwr(j)*k11(i,j)
        end do
        k11(i,i) = ki1(i) - sum
        sum = 0.0
        do j=1,nr
          sum = sum + rwr(j)*k11w(i,j)
        end do
        k11w(i,i) = ki1w(i) - sum
        if (.not. onespc) then
          sum = 0.0
          do j=1,nr
            sum = sum + rwr(j)*k22(i,j)
          end do
          k22(i,i) = ki2(i) - sum
          sum = 0.0
          do j=1,nr
            sum = sum + rwr(j)*k22w(i,j)
          end do
          k22w(i,i) = ki2w(i) - sum
        end if
      end do
c
c          Now complete the kernels with the factors rho, b1v and
c          the weights.  (The diagonal terms need only the factor b1v).
c          We use K(p = r(j),r=r(i)) = K(j,i)
c
      do   i=1,nr
        do j=1,nr
          fact = b1v(j)
          if (j .ne. i)  fact = fact*rwr(j)
          k11 (j,i) = k11 (j,i)*fact
          k11w(j,i) = k11w(j,i)*fact
          if (.not. onespc) then
            fact = b2v(j)
            if (j .ne. i)  fact = fact*rwr(j)
            k22 (j,i) = k22 (j,i)*fact
            k22w(j,i) = k22w(j,i)*fact
          end if
        end do
      end do
c
c ----------------------------------------------------------------------
c
c     Compute partial solutions, energies and fluxes
c
c     Compute the partial densities due to the sources:
c     Surface sources produce h1 and h2
c     Volume sources produce phi1 and phi2
c
      call nusrc(vthw1,swall1,vthw2,swall2,vn1,vn2)
c
c     Solve for the three partial solutions:
c     pn11,pn21 - source h1 only
c     pn12,pn22 - source h2 only
c     pn1v,pn2v - volume sources phi1 and phi2 only
c     Note: first digit indicates species, second character indicates
c     which partial solution.
c
      if (onespc)  go to 240
c
      do i=1,nr
        temp(i) = 0.0
      end do
c
      call nusolv(h1,temp,pn11,pn21)
      call nusolv(temp,h2,pn12,pn22)
      call nusolv(phi1,phi2,pn1v,pn2v)
      go to 250
c
c     Solve for partial solutions in one species case
c
  240 call nuslv1(h1,pn11)
      call nuslv1(phi1,pn1v)
c
c     Compute temperature (2/3 average energy) for each of the
c     partial solutions  (keV).  There is a surface source part
c     (the initializing value here), and an interior (scattered)
c     part at each point.
c
  250 if (onespc)  go to 280
c
      do i=1,nr
        pw11(i) = twall*(h1(i) + 2.0*h1w(i))
        pw21(i) = 0.0
        pw12(i) = 0.0
        pw22(i) = twall*(h2(i) + 2.0*h2w(i))
        pw1v(i) = 0.0
        pw2v(i) = 0.0
c
        do j=1,nr
          tk1 = tii(j)*(k11(j,i) + 2.0*k11w(j,i))
          tk2 = tii(j)*(k22(j,i) + 2.0*k22w(j,i))
          pw11(i) = pw11(i) + tk1 * (pn11(j) + fb12(j)*pn21(j))
          pw21(i) = pw21(i) + tk2 * (pn21(j) + fb21(j)*pn11(j))
          pw12(i) = pw12(i) + tk1 * (pn12(j) + fb12(j)*pn22(j))
          pw22(i) = pw22(i) + tk2 * (pn22(j) + fb21(j)*pn12(j))
          pw1v(i) = pw1v(i) + tk1 * (pn1v(j) + vn1(j) + fb12(j)*pn2v(j))
          pw2v(i) = pw2v(i) + tk2 * (pn2v(j) + vn2(j) + fb21(j)*pn1v(j))
        end do
c
        call nudiv (pw11(i),3.0*pn11(i))
        call nudiv (pw21(i),3.0*pn21(i))
        call nudiv (pw12(i),3.0*pn12(i))
        call nudiv (pw22(i),3.0*pn22(i))
        call nudiv (pw1v(i),3.0*pn1v(i))
        call nudiv (pw2v(i),3.0*pn2v(i))
      end do
c
      go to 310
c
c     Compute temperature for one species case
c
  280 do i=1,nr
        pw11(i) = twall*(h1(i) + 2.0*h1w(i))
        pw1v(i) = 0.0
        do j=1,nr
          tk1 = tii(j)*(k11(j,i) + 2.0*k11w(j,i))
          pw11(i) = pw11(i) + tk1*pn11(j)
          pw1v(i) = pw1v(i) + tk1 * (pn1v(j) + vn1(j))
        end do
        call nudiv (pw11(i),3.0*pn11(i))
        call nudiv (pw1v(i),3.0*pn1v(i))
      end do
c
c     Compute the outgoing fluxes at the surface for each of
c     the partial solutions  (1/cm**2-sec).  There is a source
c     term (shine-through of unscattered neutrals) and a
c     contribution from each of the interior points.
c     NUFLUX sets up the integral appropriate for fluxes.
c
  310 call nuflux(vthw1,vthw2)
      if (onespc)  go to 330
      pf11 = flux1*flxh1
      pf21 = 0.0
      pf12 = 0.0
      pf22 = flux2i*flxh2
      pf1v = 0.0
      pf2v = 0.0
c
      do j=1,nr
        pf11 = pf11 + flxi1(j) * (pn11(j) + fb12(j)*pn21(j))
        pf21 = pf21 + flxi2(j) * (pn21(j) + fb21(j)*pn11(j))
        pf12 = pf12 + flxi1(j) * (pn12(j) + fb12(j)*pn22(j))
        pf22 = pf22 + flxi2(j) * (pn22(j) + fb21(j)*pn12(j))
        pf1v = pf1v + flxi1(j) * (pn1v(j) + vn1(j) + fb12(j)*pn2v(j))
        pf2v = pf2v + flxi2(j) * (pn2v(j) + vn2(j) + fb21(j)*pn1v(j))
      end do
c
      go to 350
c
c     Compute fluxes for one species case
c
  330 pf11 = flux1*flxh1
      pf1v = 0.0
c
      do j=1,nr
        pf11 = pf11 + flxi1(j)*pn11(j)
        pf1v = pf1v + flxi1(j) * (pn1v(j) + vn1(j))
      end do
c
c ----------------------------------------------------------------------
c
c     Compute final densities and temperatures
c
c     We impose the reflection boundary condition (rd1,rd2<>0)
c     by solving for the proper super-position of the partial
c     solutions. (See write-up for equations.)
c     First compute normalized fluxes.
c
  350 if (onespc)  go to 420
      c11 = pf11/flux1
      c12 = pf12/flux1
      c1v = pf1v/flux1
      c21 = pf21/flux2i
      c22 = pf22/flux2i
      c2v = pf2v/flux2i
c
c     Now compute the coefficients for the super-position
c
      det   = (1.0 - rd1*c11)*(1.0-rd2*c22) - rd1*c12*rd2*c21
      alph1 = (1.0 - rd2*c22)/det
      alph2 = s2*rd1*c12/det
      alphv = ((1.0-rd2*c22)*rd1*c1v + rd1*c12*rd2*c2v)/det
      beta1 = s2*(1.0-rd1*c11)/det
      beta2 = rd2*c21/det
      betav = ((1.0-rd1*c11)*rd2*c2v + rd2*c21*rd1*c1v)/det
c
c     Compute final densities
c
      do i=1,nr
        fn11(i) = alph1*pn11(i) + beta2*pn12(i)
        fn21(i) = alph1*pn21(i) + beta2*pn22(i)
        fn12(i) = alph2*pn11(i) + beta1*pn12(i)
        fn22(i) = alph2*pn21(i) + beta1*pn22(i)
        fn1v(i) = alphv*pn11(i) + betav*pn12(i) + pn1v(i)
        fn2v(i) = alphv*pn21(i) + betav*pn22(i) + pn2v(i)
      end do
c
c     Compute final temperatures
c
      do i=1,nr
        fw11(i) = (alph1*pn11(i)*pw11(i) + beta2*pn12(i)*pw12(i))
        call nudiv (fw11(i),fn11(i))
        fw21(i) = (alph1*pn21(i)*pw21(i) + beta2*pn22(i)*pw22(i))
        call nudiv (fw21(i),fn21(i))
        fw12(i) = (alph2*pn11(i)*pw11(i) + beta1*pn12(i)*pw12(i))
        call nudiv (fw12(i),fn12(i))
        fw22(i) = (alph2*pn21(i)*pw21(i) + beta1*pn22(i)*pw22(i))
        call nudiv (fw22(i),fn22(i))
        fw1v(i) = (alphv*pn11(i)*pw11(i) + betav*pn12(i)*pw12(i) +
     .             pn1v(i)*pw1v(i))
        call nudiv (fw1v(i),fn1v(i))
        fw2v(i) = (alphv*pn21(i)*pw21(i) + betav*pn22(i)*pw22(i) +
     .             pn2v(i)*pw2v(i))
        call nudiv (fw2v(i),fn2v(i))
      end do
c
      go to 440
c
c     Compute final density and temperature in one species case
c
  420 c11     = pf11/flux1
      c1v     = pf1v/flux1
      det     = 1.0 - rd1*c11
      alph1   = 1.0 / det
      alphv   = rd1*c1v/det
c
      do i=1,nr
        fn11(i) = alph1*pn11(i)
        fn1v(i) = alphv*pn11(i) + pn1v(i)
        fw11(i) = pw11(i)
        fw1v(i) = (alphv*pn11(i)*pw11(i) + pn1v(i)*pw1v(i))
        call nudiv (fw1v(i),fn1v(i))
      end do
c
c     Interpolate back to external grid
c
  440 call nuint (nr,r,fn11,nx,rext,dn11)
      call nuint (nr,r,fn1v,nx,rext,dn1v)
      if (onespc)  go to 450
      call nuint (nr,r,fn12,nx,rext,dn12)
      call nuint (nr,r,fn21,nx,rext,dn21)
      call nuint (nr,r,fn22,nx,rext,dn22)
      call nuint (nr,r,fn2v,nx,rext,dn2v)
  450 call nuint (nr,r,fw11,nx,rext,wn11)
      call nuint (nr,r,fw1v,nx,rext,wn1v)
      if (onespc)  go to 460
      call nuint (nr,r,fw12,nx,rext,wn12)
      call nuint (nr,r,fw21,nx,rext,wn21)
      call nuint (nr,r,fw22,nx,rext,wn22)
      call nuint (nr,r,fw2v,nx,rext,wn2v)
  460 return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine aspace (rtemp, atemp)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     ASPACE converts the grid at which the A1 and A2 rates are
c     specified into an equally spaced grid.  This allows the
c     the integration routines to run more easily and faster.
c     Note that this grid is kept separately from others in the problem.
c
      include '../inc/neucg.m'
c
      dimension       rtemp(*), atemp(*)
      logical         twospc, onespc
      common /nuspcs/ twospc, onespc
      common /nuatbl/ odelr, ra(nxmax), a1(nxmax), a2(nxmax), nra, nram1
c
      dr = ra(nra)/(nra-1)
      do i=1,nra
        rtemp(i) = ra(i)
        ra   (i) = dr * (i-1)
        atemp(i) = a1(i)
      end do
c
      call nuint (nra, rtemp, atemp, nra, ra, a1)
      if (onespc)  return
c
      do i=1,nra
        atemp(i) = a2(i)
      end do
      call nuint (nra, rtemp, atemp, nra, ra, a2)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nuerr (ierr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     NUERR handles the I/O for errors that are detected by NEUCG2.
c     We use standard I/O files for the ONETWO transport code
c
      data nout/7/, nqik/8/, ncrt/6/
c
      if (ierr .eq. 1)  write (nout, 9001)
      if (ierr .eq. 1)  write (nqik, 9001)
      if (ierr .eq. 1)  write (ncrt, 9001)
 9001 format (' Neucg2 ERROR: dimension nxmax too small - recompile')
      if (ierr .eq. 2)  write (nout, 9002)
      if (ierr .eq. 2)  write (nqik, 9002)
      if (ierr .eq. 2)  write (ncrt, 9002)
 9002 format (' Neucg2 ERROR: solution failed to converge')
      if (ierr .eq. 3)  write (nout, 9003)
      if (ierr .eq. 3)  write (nqik, 9003)
      if (ierr .eq. 3)  write (ncrt, 9003)
 9003 format (' Neucg2 ERROR: error in IMSL spline routines')
      if (ierr .eq. 4)  write (nout, 9004)
      if (ierr .eq. 4)  write (nqik, 9004)
      if (ierr .eq. 4)  write (ncrt, 9004)
 9004 format (' Neucg2 WHOOPS: internal error')
      if (ierr .eq. 5)  write (nout, 9005)
      if (ierr .eq. 5)  write (nqik, 9005)
      if (ierr .eq. 5)  write (ncrt, 9005)
 9005 format (' Neucg2 ERROR: flux1 cannot be zero')
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nugrid (rminor)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     NUGRID sets up the base points and weights for the radial
c     and theta quadrature:
c     radial - 22 point Lobatto, 0 to rminor.
c     theta  - 10 point Gaussian in EXP (theta), 0 to pi/2
c     the exp moves points closer to pi/2, which is singular at p = r.
c     The weights include the normalizing factor for interval length.
c
      include '../inc/neucg.m'
      dimension       vi(11), wi(11), v(5), w(5)
      common /nuquad/ r(nr), wr(nr), rwr(nr), sint(nt), wt(nt)
      data            prevrm /0.0/
c
c     these abscissae and weights are for 22-point lobatto quadrature
c
      data vi/1.0           , 0.984152438   , 0.947204284   ,
     .        0.890062290   , 0.813948928   , 0.720487240   ,
     .        0.611669438   , 0.489814875   , 0.357520710   ,
     .        0.217606585   , 0.7305454e-1/ ,
     .     wi/0.432900433e-2, 0.265457477e-1, 0.472144653e-1,
     .        0.668656059e-1, 0.850900604e-1, 0.101500575   ,
     .        0.115747645   , 0.127527697   , 0.136589689   ,
     .        0.142740492   , 0.145849019 /
c
c     these abscissae and weights are for the angular quadrature
c
      data v/0.973906529   , 0.865063367, 0.679409568, 0.433395394,
     .       0.148874339/  ,
     .     w/0.666713443e-1, 0.149451349, 0.219086363, 0.269266719,
     .       0.295524225/
c
c     check if grid is same as previous grid
c
      if (rminor .eq. prevrm)  return
c
c     set up radial quadrature
c
      rmin2 = rminor / 2.0
      k     = nr / 2
c
c     work from both ends in to middle
c
      do i=1,k
        r(i)   = rmin2*(1.0-vi(i))
        wr(i)  = rmin2*wi(i)
        rwr(i) = r(i)*wr(i)
        j      = nr + 1 - i
        r(j)   = rmin2 * (1.0 + vi(i))
        wr(j)  = wr(i)
        rwr(j) = r(j) * wr(j)
      end do
c
      if (prevrm .ne. 0.0)  go to 120
c
c     set up theta quadrature
c
      k = nt / 2
c
c     work from both ends in to middle
c
      do i=1,k
        theta   = LOG (1.90523869*(-v(i)) + 2.90523869)
        sint(i) = SIN (theta)
        wt(i)   = w(i)/(-v(i)+1.52486862)
        j       = nt + 1 - i
        theta   = LOG (1.90523869*v(i) + 2.90523869)
        sint(j) = SIN (theta)
        wt(j)   = w(i) / (v(i) + 1.52486862)
      end do
c
  120 prevrm = rminor
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nuint (n1, r1, f1, n2, r2, f2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          NUINT converts the array of values f1, specified at radii r1,
c          into values f2, corresponding to the given radii r2.
c          To prevent getting negative values, the fitting is done
c          on the log of the values f1.  The IMSL cubic spline routine
c          is used to interpolate the logarithms.
c
      include '../inc/neucg.m'
      dimension  r1(*), f1(*), r2(*), f2(*)
      dimension  bpar(4), c(nxmax,3), f1log(nxmax)
      data       bpar /4*0.0/
c
      do 100 i=1,n1
      if (f1(i) .le. 0.0)  go to 120
  100 f1log(i) = LOG (f1(i))
c
      call icsicu1 (r1, f1log, n1, bpar, c, nxmax, ier)
      if (ier .ne. 0)  call nuerr(3)
      call icsevu1 (r1, f1log, n1, c, nxmax, r2, f2, n2, ier)
      if (ier .ne. 0)  call nuerr(3)
c
      do 110 i=1,n2
  110 f2(i) = EXP (f2(i))
      return
c
c     Some elements are non-positive: use standard spline
c
  120 continue
      call icsicu1 (r1, f1, n1, bpar, c, nxmax, ier)
      if (ier .ne. 0)  call nuerr(3)
      call icsevu1 (r1, f1, n1, c, nxmax, r2, f2, n2, ier)
      if (ier .ne. 0)  call nuerr(3)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nusrc (vthw1, swall1, vthw2, swall2, vn1, vn2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          NUSRC calculates the source density functions:
c          h1 and h2 for the surface source, and phi1 and phi2 for
c          the volume source.  H1w and h2w are used later in the energy
c          computation.  These source densities are then the starting
c          guess in the iterative solution of the integral equations.
c
      include '../inc/neucg.m'
      dimension       vn1(*), vn2(*)
      logical         twospc, onespc
      common /nuspcs/ twospc, onespc
      common /nuquad/ r(nr),wr(nr),rwr(nr),sint(nt),wt(nt)
      common /nuatr0/ atrr01(nr,nt),atar01(nr,nt),
     .                atrr02(nr,nt),atar02(nr,nt)
      common /nukern/ fb12(nr),fb21(nr),k11(nr,nr),k11w(nr,nr),
     .                                  k22(nr,nr),k22w(nr,nr)
      real*8 k11,k11w,k22,k22w
      common /nuhphi/ h1(nr),h1w(nr),phi1(nr),h2(nr),h2w(nr),phi2(nr)
c
      ovw1 = 1.0 / vthw1
      if (twospc)  ovw2 = 1.0 / vthw2
c
c          Compute surface source density
c
      do 110 i=1,nr
      sum1 = 0.0
      sum1w = 0.0
      sum2 = 0.0
      sum2w = 0.0
      do 100 k=1,nt
c
c          do first species
c          A+ and A- are found from previously computed At's
c
      atr = atrr01(i,k)
      ata = atar01(i,k)
      apovw = (ata + atr)*ovw1
      amovw = (ata - atr)*ovw1
      g1ap = g1(apovw)
      g1am = g1(amovw)
      g3ap = g1ap + 0.5*apovw*g0(apovw)
      g3am = g1am + 0.5*amovw*g0(amovw)
      sum1 = sum1 + wt(k)*(g1ap + g1am)
      sum1w = sum1w + wt(k)*(g3ap + g3am)
      if (onespc)  go to 100
c
c          do second species
c
      atr = atrr02(i,k)
      ata = atar02(i,k)
      apovw = (ata + atr)*ovw2
      amovw = (ata - atr)*ovw2
      g1ap = g1(apovw)
      g1am = g1(amovw)
      g3ap = g1ap + 0.5*apovw*g0(apovw)
      g3am = g1am + 0.5*amovw*g0(amovw)
      sum2 = sum2 + wt(k)*(g1ap + g1am)
      sum2w = sum2w + wt(k)*(g3ap + g3am)
  100 continue
      h1(i) = swall1*sum1
      h1w(i) = swall1*sum1w
      if (onespc)  go to 110
      h2(i) = swall2*sum2
      h2w(i) = swall2*sum2w
  110 continue
c
c          Now compute the volume source density.
c          Note that at this point the kernels are complete
c
c  The following compiler directives were added to allow object code
c     compiled with cft 1.11 on a CRAY-1 to execute properly on a
c     CRAY X-MP.
c
      do 200 i=1,nr
      phi1(i) = 0.0
      do 200 j=1,nr
      phi1(i) = phi1(i) + k11(j,i)*vn1(j)
  200 continue
      if (onespc)  go to 220
c
      do 210 i=1,nr
      phi2(i) = 0.0
      do 210 j=1,nr
      phi2(i) = phi2(i) + k22(j,i)*vn2(j)
  210 continue
c
  220 continue
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nusolv (sn1, sn2, pn1, pn2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          NUSOLV solves the 2 species coupled integral equations.
c          The technique is Gauss-Seidel iteration.
c          Note that the diagonal terms of K(i,j) are non-zero now.
c          To get the equations into the standard form, we subtract
c          the K(i,i) term from both sides and then divide by 1.0-K(i,i).
c
      include '../inc/neucg.m'
      dimension       sn1(*),sn2(*),pn1(*),pn2(*)
      dimension       f1(nr),f2(nr)
      common /nukern/ fb12(nr),fb21(nr),k11(nr,nr),k11w(nr,nr),
     .                                  k22(nr,nr),k22w(nr,nr)
      real*8          k11, k11w, k22, k22w
      data            maxit/50/, tol/1.0e-4/
c
c          Set up factor with diagonal term, and initial guess
c
      do 100 i=1,nr
      f1(i) = 1.0 / (1.0-k11(i,i))
      f2(i) = 1.0 / (1.0-k22(i,i))
      pn1(i) = sn1(i)
      pn2(i) = sn2(i)
  100 continue
c
c          Now iterate until the solution converges
c
      do 150 it=1,maxit
      del = 0.0
      do 120 i=1,nr
      tn1 = sn1(i)
      do 110 j=1,nr
      if (j .eq. i)  go to 110
      tn1 = tn1 + k11(j,i)*(pn1(j) + fb12(j)*pn2(j))
  110 continue
      tn1 = tn1 + k11(i,i)*fb12(i)*pn2(i)
      tn1 = tn1*f1(i)
      del1 = 0.0
      if (pn1(i) .ne. 0.0) del1 = ABS ((tn1-pn1(i))/pn1(i))
      if (pn1(i) .eq. 0.0 .and. tn1 .ne. 0.0) del1 = 1.0
      del = MAX (del, del1)
      pn1(i) = tn1
  120 continue
c
      do 140 i=1,nr
      tn2 = sn2(i)
      do 130 j=1,nr
      if (j .eq. i)  go to 130
      tn2 = tn2 + k22(j,i)*(pn2(j) + fb21(j)*pn1(j))
  130 continue
      tn2 = tn2 + k22(i,i)*fb21(i)*pn1(i)
      tn2 = tn2*f2(i)
      del2 = 0.0
      if (pn2(i) .ne. 0.0) del2 = ABS ((tn2-pn2(i))/pn2(i))
      if (pn2(i) .eq. 0.0 .and. tn2 .ne. 0.0) del2 = 1.0
      del = MAX (del,del2)
      pn2(i) = tn2
  140 continue
      if (del .lt. tol)  return
  150 continue
      call nuerr(2)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nuslv1 (sn1, pn1)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          NUSLV1 solves the 1 species integral equation.
c          The technique is Gauss-Seidel iteration.
c          Note that the diagonal terms of K(i,j) are non-zero now.
c          To get the equations into the standard form, we subtract
c          the K(i,i) term from both sides and then divide by 1.0-K(i,i).
c
      include '../inc/neucg.m'
      dimension       sn1(*), pn1(*)
      dimension       f1(nr)
      common /nukern/ fb12(nr),fb21(nr),k11(nr,nr),k11w(nr,nr),
     .                                  k22(nr,nr),k22w(nr,nr)
      real*8          k11,k11w,k22,k22w
      data            maxit/50/,tol/1.0e-4/
c
c     Set up factor with diagonal term, and initial guess
c
      do i=1,nr
        f1 (i) = 1.0 / (1.0-k11(i,i))
        pn1(i) = sn1(i)
      end do
c
c     Now iterate until the solution converges
c
      do 130 it=1,maxit
      del = 0.0
      do 120 i=1,nr
      tn1 = sn1(i)
      do 110 j=1,nr
      if (j .eq. i)  go to 110
      tn1 = tn1 + k11(j,i)*pn1(j)
  110 continue
      tn1 = tn1*f1(i)
      del1 = 0.0
      if (pn1(i) .ne. 0.0) del1 = ABS ((tn1-pn1(i))/pn1(i))
      if (pn1(i) .eq. 0.0 .and. tn1 .ne. 0.0) del1 = 1.0
      del    = MAX (del, del1)
      pn1(i) = tn1
  120 continue
      if (del .lt. tol)  return
  130 continue
      call nuerr(2)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nudiv (w, dn)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     NUDIV prevents a division by zero if a density is zero
c
      if (dn .eq. 0.0)  w = 0.0
      if (dn .ne. 0.0)  w = w / dn
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine nuflux (vthw1, vthw2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          NUFLUX sets up the integrals needed to compute the outward
c          going neutral flux at the surface for the partial solutions.
c          There is a source term, which can be thought of as the
c          shine-through of unscattered neutrals, and a contribution
c          from each of the interior points.
c          Note that we just set up the integrals here based on the
c          inputs - there is no use of information from the solutions.
c
      include '../inc/neucg.m'
      logical         twospc, onespc
      common /nuspcs/ twospc, onespc
      common /nuquad/ r(nr),wr(nr),rwr(nr),sint(nt),wt(nt)
      common /nuatr0/ atrr01(nr,nt),atar01(nr,nt),
     .                atrr02(nr,nt),atar02(nr,nt)
      common /nuovth/ ovth1(nr),ovth2(nr)
      common /nubvth/ b1v(nr),b2v(nr),vth1(nr),vth2(nr)
      common /nuflxi/ flxh1,flxh2,flxi1(nr),flxi2(nr)
c
      cost(k) = SQRT (1.0-sint(k)**2)
c
c          Compute surface source term
c
      rminor = r(nr)
      flxh1 = 0.0
      if (twospc) flxh2 = 0.0
      do 100 k=1,nt
      flxh1 = flxh1 + wt(k)*cost(k)*g2(2.0*atar01(nr,k)/vthw1)
  100 continue
      flxh1 = flxh1*2.
      if (onespc)  go to 120
      do 110 k=1,nt
      flxh2 = flxh2 + wt(k)*cost(k)*g2(2.0*atar02(nr,k)/vthw2)
  110 continue
      flxh2 = flxh2*2.
  120 continue
c
c     Compute contribution from each interior point
c
      do 140 j=1,nr
      sum1 = 0.0
      sum2 = 0.0
      do 130 k=1,nt
c
c          do first species
c
      atp = atrr01(j,k)
      ata = atar01(j,k)
      apovp = (ata + atp)*ovth1(j)
      amovp = (ata - atp)*ovth1(j)
      g1ap = g1(apovp)
      g1am = g1(amovp)
      sum1 = sum1 + wt(k)*(g1ap + g1am)
      if (onespc)  go to 130
c
c          do second species
c
      atp = atrr02(j,k)
      ata = atar02(j,k)
      apovp = (ata + atp)*ovth2(j)
      amovp = (ata - atp)*ovth2(j)
      g1ap = g1(apovp)
      g1am = g1(amovp)
      sum2 = sum2 + wt(k)*(g1ap + g1am)
  130 continue
      flxi1(j) = rwr(j)*b1v(j)*vth1(j)*sum1/rminor
      if (onespc)  go to 140
      flxi2(j) = rwr(j)*b2v(j)*vth2(j)*sum2/rminor
  140 continue
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine kernel
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     KERNEL sets up the kernels for the one- or two-species case.
c     Only the theta integral is calculated, the factors b1v(p) etc.
c     are multiplied in later. Also note that the factor rho is
c     left out until the rho integral is actually performed.
c
      include '../inc/neucg.m'
      logical         twospc, onespc
      common /nuspcs/ twospc, onespc
      common /nuquad/ r(nr),wr(nr),rwr(nr),sint(nt),wt(nt)
      common /nuatr0/ atrr01(nr,nt),atar01(nr,nt),
     .                atrr02(nr,nt),atar02(nr,nt)
      common /nuovth/ ovth1(nr),ovth2(nr)
      common /nukern/ fb12(nr),fb21(nr),k11(nr,nr),k11w(nr,nr),
     .                                  k22(nr,nr),k22w(nr,nr)
      real*8          k11, k11w, k22, k22w
      data            pi /3.141592654/
c
c          Set up kernels. Note that w denotes the energy kernel.
c          The second subscript of K determines the fixed radius r at
c          which the rho integral will be evaluated.  The first subscript
c          determines the base point rho of the quadrature (at r) - it
c          corresponds to the point used for vth(rho).  Hence we have
c          K(p = r(i),r=r(j)) = K(i,j); and K(r=r(j),p=r(i)) = K(j,i).
c          Note that since only rholt and rhogt appear elsewhere in
c          the kernel, and since the fixed radii are the quadrature
c          points, we can compute both rholt = r, rhogt=p and
c          rholt = p, rhogt=r at the same time.
c          Thus the first loop on i is the rholt loop, and the second
c          loop on j is the rhogt loop:  rholt = r(i), rhogt=r(j).
c
      do 130 i=1,nr
      k11 (i,i) = 0.0
      k11w(i,i) = 0.0
      if (twospc) k22 (i,i) = 0.0
      if (twospc) k22w(i,i) = 0.0
c
c          Compute the At functions.  atrr01(j,k) =  At at r(j) from the
c          path starting at r0 = r(i)*sint(k).  atar01 saves the value
c          at rminor for use later.  Note that after the outer loop on i
c          is done, atrr01(i,k)  =   At at r(i) from r0=r(i)*sint(k), also
c          to be used later.
c
      do 100 k=1,nt
      r0 = r(i)*sint(k)
      if (twospc)  call atv2(r(i),nr-i+1,r0,atrr01(i,k),atrr02(i,k))
      if (onespc)  call atv1(r(i),nr-i+1,r0,atrr01(i,k))
      atar01(i,k) = atrr01(nr,k)
      if (twospc) atar02(i,k) = atrr02(nr,k)
      if (i .ne. 1)  go to 100
c
c          When i = 1, rholt=0, and there is no dependence on theta(k).
c          So we can copy the k = 1 values into all the k=2,nt slots.
c
      do kk=2,nt
        do ii=1,nr
          atrr01(ii,kk) = atrr01(ii,1)
          if (twospc) atrr02(ii,kk) = atrr02(ii,1)
        end do
        atar01(1,kk) = atrr01(nr,kk)
        if (twospc) atar02(1,kk) = atrr02(nr,kk)
      end do
c
      go to 105
  100 continue
c
  105 if (i .eq. nr)  go to 130
      j0 = i+1
      do 120 j=j0,nr
      if (i .eq. 1)  go to 111
      sum1p  = 0.0
      sum1r  = 0.0
      sum1wp = 0.0
      sum1wr = 0.0
      sum2p  = 0.0
      sum2r  = 0.0
      sum2wp = 0.0
      sum2wr = 0.0
      pgt2   = r(j) * r(j)
      do 110 k=1,nt
      r0     = r(i)*sint(k)
      wtx    = wt(k) / SQRT (pgt2 - r0*r0)
c
c          do first species
c          here ap = A+,  am = A-
c
      atplt = atrr01(i,k)
      atpgt = atrr01(j,k)
      ap    = atpgt + atplt
      am    = atpgt - atplt
c
c          first: r = rholt=r(i), p=rhogt=r(j)
c
      apovp = ap*ovth1(j)
      amovp = am*ovth1(j)
      call g0g2(apovp,g0ap,g2ap)
      call g0g2(amovp,g0am,g2am)
      sum1p = sum1p + wtx*(g0ap + g0am)
      sum1wp = sum1wp + wtx*(g2ap + g2am)
c
c          now: r = rhogt=r(j), p=rholt=r(i)
c
      apovp = ap*ovth1(i)
      amovp = am*ovth1(i)
      call g0g2(apovp,g0ap,g2ap)
      call g0g2(amovp,g0am,g2am)
      sum1r = sum1r + wtx*(g0ap + g0am)
      sum1wr = sum1wr + wtx*(g2ap + g2am)
      if (onespc)  go to 110
c
c          do second species
c
      atplt = atrr02(i,k)
      atpgt = atrr02(j,k)
      ap = atpgt + atplt
      am = atpgt - atplt
c
c          first: r = rholt=r(i), p=rhogt=r(j)
c
      apovp = ap*ovth2(j)
      amovp = am*ovth2(j)
      call g0g2(apovp,g0ap,g2ap)
      call g0g2(amovp,g0am,g2am)
      sum2p = sum2p + wtx*(g0ap + g0am)
      sum2wp = sum2wp + wtx*(g2ap + g2am)
c
c          now: r = rhogt=r(j), p=rholt=r(i)
c
      apovp = ap*ovth2(i)
      amovp = am*ovth2(i)
      call g0g2(apovp,g0ap,g2ap)
      call g0g2(amovp,g0am,g2am)
      sum2r = sum2r + wtx*(g0ap + g0am)
      sum2wr = sum2wr + wtx*(g2ap + g2am)
  110 continue
      k11(j,i) = sum1p
      k11(i,j) = sum1r
      k11w(j,i) = sum1wp
      k11w(i,j) = sum1wr
      if (onespc)  go to 120
      k22(j,i) = sum2p
      k22(i,j) = sum2r
      k22w(j,i) = sum2wp
      k22w(i,j) = sum2wr
      go to 120
c
c          When i = 1, rholt=0, and there is no dependence on theta(k).
c          So we can do that case very quickly here.
c
  111 apm = atrr01(j,1)
      apovp = apm*ovth1(j)
      call g0g2(apovp,g0ap,g2ap)
      k11(j,1) = pi*g0ap/r(j)
      k11w(j,1) = pi*g2ap/r(j)
      apovp = apm*ovth1(1)
      call g0g2(apovp,g0ap,g2ap)
      k11(1,j) = pi*g0ap/r(j)
      k11w(1,j) = pi*g2ap/r(j)
      if (onespc)  go to 120
c
      apm = atrr02(j,1)
      apovp = apm*ovth2(j)
      call g0g2(apovp,g0ap,g2ap)
      k22(j,1) = pi*g0ap/r(j)
      k22w(j,1) = pi*g2ap/r(j)
      apovp = apm*ovth2(1)
      call g0g2(apovp,g0ap,g2ap)
      k22(1,j) = pi*g0ap/r(j)
      k22w(1,j) = pi*g2ap/r(j)
  120 continue
  130 continue
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine kintegral
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c    (KINTEGRAL used to be named KINT but was changed on 28 March 1999
c     by Joe Freeman to avoid a name conflict with the KINT intrinsic.)
c
c     KINTEGRAL computes the singular integral using the quadrature
c     method due to Krylov and Paltsev.
c     The integral at each point r is over rho times K(r,p).
c     The factors b1v(r) etc. are multiplied in later.
c     Note that at r = 0 the integral is not singular,
c     and this method does not apply.
c
      include '../inc/neucg.m'
      parameter      (nkint = 10)
      logical         twospc, onespc
      common /nuspcs/ twospc, onespc
      common /nuquad/ r(nr),wr(nr),rwr(nr),sint(nt),wt(nt)
      common /nukint/ ki1(nr),ki1w(nr),ki2(nr),ki2w(nr)
      real*8          ki1,    ki1w,    ki2,    ki2w
      dimension       rk(nkint), wk(nkint)
      data            pi /3.141592654/
c
c     r = 0 is a special case
c
      ki1 (1) = wr(1) * pi
      ki1w(1) = wr(1) * pi / 2.0
      if (twospc)  ki2 (1) = ki1 (1)
      if (twospc)  ki2w(1) = ki1w(1)
c
c     compute integral at each point r = r(i)
c
      do 130 i=2,nr
      sum1  = 0.0
      sum1w = 0.0
      sum2  = 0.0
      sum2w = 0.0
      rr    = r(i)
      amr   = r(nr) - rr
c
c     do 0 to r(i) segment
c
      call kgrid(1,rr,amr,rk,wk)
c
c     note that wk includes a factor rho
c     and that rk(ik) = rholt < rr = rhogt
c
      do 100 ik=1,nkint
        call kintth(rk(ik),rr,i,ti1,ti1w,ti2,ti2w)
        sum1  = sum1 + wk(ik)*ti1
        sum1w = sum1w + wk(ik)*ti1w
        if (onespc)  go to 100
        sum2  = sum2 + wk(ik)*ti2
        sum2w = sum2w + wk(ik)*ti2w
  100 continue
c
      ki1(i)  = rr*sum1
      ki1w(i) = rr*sum1w
      if (onespc)  go to 110
      ki2(i)  = rr*sum2
      ki2w(i) = rr*sum2w
c
c     now do r(i) to rminor seqment
c
  110 if (i .eq. nr)  go to 130
      sum1  = 0.0
      sum1w = 0.0
      sum2  = 0.0
      sum2w = 0.0
      call kgrid(2,rr,amr,rk,wk)
c
c     note that rr = rholt < rk(ik) = rhogt
c
      do 120 ik=1,nkint
        call kintth(rr,rk(ik),i,ti1,ti1w,ti2,ti2w)
        sum1  = sum1 + wk(ik)*ti1
        sum1w = sum1w + wk(ik)*ti1w
        if (onespc)  go to 120
        sum2  = sum2 + wk(ik)*ti2
        sum2w = sum2w + wk(ik)*ti2w
  120 continue
c
      ki1(i)  = ki1(i) + amr*sum1
      ki1w(i) = ki1w(i) + amr*sum1w
      if (onespc)  go to 130
      ki2(i)  = ki2(i) + amr*sum2
      ki2w(i) = ki2w(i) + amr*sum2w
  130 continue
c
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function g0 (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      g0f1(x) =
     .(9.9995058e-1 + x*(7.7628443e+2 + x*(-4.3649686e+3 +
     . x*6.1480022e+4)))/
     .(1.0 + x*7.8621464e+2)
      g0f2(x) =
     .(9.9703418e-1 + x*(7.7783588e+1 + x*(3.9012546e+2 +
     . x*(-8.9205431e+2 + x*1.0037943e+3))))/
     .(1.0 + x*(8.4398466e+1 + x*7.1649574e+2))
      g0f3(x) =
     .(9.7552510e-1 + x*(7.0647154 + x*(-4.0953920 +
     . x*(9.0556774e-1 - x*5.8788928e-2))))/
     .(1.0 + x*(1.1344206e+1 + x*1.5956867e+1))
      g0f4(x) =
     .(8.4513992e-1 + x*(-2.2811875e-1 + x*(2.5926818e-2 +
     . x*(-1.4549910e-3 + x*3.3570582e-5))))/
     .(1.0 + x*(2.0520190  + x*(6.1145865e-1 + x*1.2572764e-1)))
      g0f5(x) =
     .(1.9937501e-3 + x*(-3.2073160e-4 + x*(1.7925104e-5 -
     . x*3.4571807e-7)))/
     .(1.0 + x*(-3.3316230e-1 + x*3.6461690e-2))
      g0f6(x) =
     .(-2.4903487e-7 + x*7.3163546e-9)/
     .(1.0 + x*(-3.2915104e-1 + x*(4.3949080e-2 + x*(-2.9908526e-3 +
     . x*(1.0457349e-4 - x*1.5316628e-6)))))
c
c          the above fits were obtained with the IMSL routine iratcu,
c          and typical max. relative error is 3.e-5.
c
      if (x .gt. 8.0)  go to 11
      if (x .lt. 1.0)  go to 14
c
c          1 < x < 8
c
      g0 = g0f4(x)
      return
   11 if (x .gt. 15.0)  go to 12
c
c          8 < x < 15
c
      g0 = g0f5(x)
      return
   12 if (x .gt. 23.0)  go to 13
c
c          15 < x < 23
c
      g0 = g0f6(x)
      return
c
c       23 < x  - cut off the tails is ok, because the g's are down
c                 typically by x1e6 or more, and are added to a much larger term
c
   13 g0 = 0.0
      return
   14 if (x .lt. 0.1)  go to 15
c
c          0.1 < x < 1
c
      g0 = g0f3(x)
      return
   15 if (x .lt. 0.01)  go to 16
c
c          0.01 < x < 0.1
c
      g0 = g0f2(x)
      return
c
c          0.0 < x < 0.01
c
   16 g0 = g0f1(x)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function g1 (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      g1f1(x) =
     .(5.6418256e-1 + x*(7.6987602 + x*(1.7382026e-1 +
     . x*(-3.7205545 + x*(2.3330170 - x*5.5200657e-1)))))/
     .(1.0 + x*(1.5405905e+1 + x*2.2881654e+1))
      g1f2(x) =
     .(5.4679544e-1 + x*(-1.3486195e-1 + x*(1.3867652e-2 +
     . x*(-6.9690004e-4 + x*1.4249435e-5))))/
     .(1.0 + x*(1.2227160 + x*(2.9519202e-1 + x*5.2677797e-2)))
      g1f3(x) =
     .(-2.8658332e-3 + x*(2.1641101e-4 - x*4.3331406e-6))/
     .(1.0 + x*(-6.7002968e-1 + x*(1.7850153e-1 + x*(-2.6814973e-2 +
     . x*(1.9753437e-3 - x*8.6337343e-5)))))
      g1f4(x) =
     .(5.2280034e-7 - x*1.3739458e-8)/
     .(1.0 + x*(-3.8248068e-1 + x*(6.1669795e-2 + x*(-5.3850087e-3 +
     . x*(2.7014959e-4 + x*(-7.4475777e-6 + x*9.0283061e-8))))))
c
c          the above fits were obtained with the IMSL routine iratcu,
c          and typical max. relative error is 3.e-5.
c
      if (x .gt. 8.0)  go to 11
      if (x .lt. 1.0)  go to 14
c
c          1.<x<8.
c
      g1 = g1f2(x)
      return
   11 if (x .gt. 15.0)  go to 12
c
c          8.<x<15.
c
      g1 = g1f3(x)
      return
   12 if (x .gt. 23.0)  go to 13
c
c          15.<x<23.
c
      g1 = g1f4(x)
      return
c          23.<x  - cut off the tails is ok, because the g's are down
c                   typically by x1e6 or more, and are added to a much
c                   larger term.
   13 g1 = 0.0
      return
c
c          0.<x<1.
c
   14 g1 = g1f1(x)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function g2 (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      g2f1(x) =
     .(4.9999949e-1 + x*(1.6514292 + x*(-7.7504797e-1 +
     . x*(1.3080537e-1 - x*4.4677504e-3))))/
     .(1.0 + x*(4.4309042 + x*2.4732245))
      g2f2(x) =
     .(4.9648578e-1 + x*(-1.1727865e-1 + x*(1.1589727e-2 +
     . x*(-5.6263257e-4 + x*1.1186908e-5))))/
     .(1.0 + x*(8.4868414e-1 + x*(1.7879136e-1 + x*2.5181820e-2)))
      g2f3(x) =
     .(3.1672895e-3 + x*(-2.6068900e-4 +x*5.6722220e-6))/
     .(1.0 + x*(-5.2547035e-1 + x*(1.1447531e-1 + x*(-1.1069281e-2 +
     . x*5.7336289e-4))))
      g2f4(x) =
     .(1.8416390e-6 - x*4.8297569e-8)/
     .(1.0 + x*(-3.9004312e-1 + x*(6.4195056e-2 + x*(-5.7286242e-3 +
     . x*(2.9425121e-4 + x*(-8.3255466e-6 + x*1.0435769e-7))))))
c
c          the above fits were obtained with the IMSL routine iratcu,
c          and typical max. relative error is 3.e-5.
c
      if (x .gt. 8.0)  go to 11
      if (x .lt. 1.0)  go to 14
c
c          1.<x<8.
c
      g2 = g2f2(x)
      return
   11 if (x .gt. 15.0)  go to 12
c
c          8.<x<15.
c
      g2 = g2f3(x)
      return
   12 if (x .gt. 23.0)  go to 13
c
c          15.<x<23.
c
      g2 = g2f4(x)
      return
c          23.<x  - cut off the tails is ok, because the g's are down
c                   typically by x1e6 or more, and are added to a much
c                   larger term.
   13 g2 = 0.0
      return
c
c          0.<x<1.
c
   14 g2 = g2f1(x)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine g0g2 (x, g0, g2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          g0g2 returns both g0 and g2, as is usually desired.
c          timing tests show this combined call runs about 1.3x
c          faster than separate calls (which is important here).
c
      g0f1(x) =
     .(9.9995058e-1 + x*(7.7628443e+2 + x*(-4.3649686e+3 +
     . x*6.1480022e+4)))/
     .(1.0 + x*7.8621464e+2)
      g0f2(x) =
     .(9.9703418e-1 + x*(7.7783588e+1 + x*(3.9012546e+2 +
     . x*(-8.9205431e+2 + x*1.0037943e+3))))/
     .(1.0 + x*(8.4398466e+1 + x*7.1649574e+2))
      g0f3(x) =
     .(9.7552510e-1 + x*(7.0647154 + x*(-4.0953920 +
     . x*(9.0556774e-1 - x*5.8788928e-2))))/
     .(1.0 + x*(1.1344206e+1 + x*1.5956867e+1))
      g0f4(x) =
     .(8.4513992e-1 + x*(-2.2811875e-1 + x*(2.5926818e-2 +
     . x*(-1.4549910e-3 + x*3.3570582e-5))))/
     .(1.0 + x*(2.0520190  + x*(6.1145865e-1 + x*1.2572764e-1)))
      g0f5(x) =
     .(1.9937501e-3 + x*(-3.2073160e-4 + x*(1.7925104e-5 -
     . x*3.4571807e-7)))/
     .(1.0 + x*(-3.3316230e-1 + x*3.6461690e-2))
      g0f6(x) =
     .(-2.4903487e-7 + x*7.3163546e-9)/
     .(1.0 + x*(-3.2915104e-1 + x*(4.3949080e-2 + x*(-2.9908526e-3 +
     . x*(1.0457349e-4 - x*1.5316628e-6)))))
      g2f1(x) =
     .(4.9999949e-1 + x*(1.6514292 + x*(-7.7504797e-1 +
     . x*(1.3080537e-1 - x*4.4677504e-3))))/
     .(1.0 + x*(4.4309042 + x*2.4732245))
      g2f2(x) =
     .(4.9648578e-1 + x*(-1.1727865e-1 + x*(1.1589727e-2 +
     . x*(-5.6263257e-4 + x*1.1186908e-5))))/
     .(1.0 + x*(8.4868414e-1 + x*(1.7879136e-1 + x*2.5181820e-2)))
      g2f3(x) =
     .(3.1672895e-3 + x*(-2.6068900e-4 +x*5.6722220e-6))/
     .(1.0 + x*(-5.2547035e-1 + x*(1.1447531e-1 + x*(-1.1069281e-2 +
     . x*5.7336289e-4))))
      g2f4(x) =
     .(1.8416390e-6 - x*4.8297569e-8)/
     .(1.0 + x*(-3.9004312e-1 + x*(6.4195056e-2 + x*(-5.7286242e-3 +
     . x*(2.9425121e-4 + x*(-8.3255466e-6 + x*1.0435769e-7))))))
c
c          the above fits were obtained with the IMSL routine iratcu,
c          and typical max. relative error is 3.e-5.
c
      if (x .gt. 8.0)  go to 11
      if (x .lt. 1.0)  go to 14
c
c          1 < x < 8
c
      g0 = g0f4(x)
      g2 = g2f2(x)
      return
   11 if (x .gt. 15.0)  go to 12
c
c          8 < x < 15
c
      g0 = g0f5(x)
      g2 = g2f3(x)
      return
   12 if (x .gt. 23.0)  go to 13
c
c          15 < x < 23
c
      g0 = g0f6(x)
      g2 = g2f4(x)
      return
c          23 < x - cut off the tails is ok, because the g's are down
c                   typically by x1e6 or more, and are added to a much
c                   larger term.
   13 g0 = 0.0
      g2 = 0.0
      return
   14 if (x .lt. 0.1)  go to 15
c
c          0.1 < x < 1
c
      g0 = g0f3(x)
      g2 = g2f1(x)
      return
   15 if (x .lt. 0.01)  go to 16
c
c          0.01 < x < 0.1
c
      g0 = g0f2(x)
      g2 = g2f1(x)
      return
c
c          0 < x < 0.01
c
   16 g0 = g0f1(x)
      g2 = g2f1(x)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine atv1 (r, n, r0, at1)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     Atv1 calculates the value of the At integral from
c     r0 = rholt*SIN (theta) out to each of the nr radii specified
c     in r. This saves on redundant integration along the same
c     path. The integration is done by the trapezoidal rule.
c
      include '../inc/neucg.m'
      dimension       r(*), at1(*)
      real*8          singloid ! single precision variable
      common /nuatbl/ odelr, ra(nxmax), a1(nxmax), a2(nxmax), nra, nram1
c
c     This statement function interpolates the A1 rate
c     to the desired point.  It assumes an equally spaced grid with
c     ra(1) = 0 and ra(2) = delta r, as produced by the routine Aspace.
c     (odelr = 1/delta r)
c     Note that x is a leg of the right triangle with hypotenuse
c     rho and the other leg r0, thus following the neutral's path.
c
      a1f(i,rdum) = a1(i) + (a1(i+1)-a1(i)) * (rdum-ra(i))*odelr
      x(rdummy)   = SQRT (rdummy*rdummy - r0sq)
c
c     We start at r0
c
      r0sq     = r0 * r0
      a1int    = 0.0
      singloid = r0 * odelr
      i        = MIN0 (IFIX (singloid) + 1, nram1)
c
c     aij: i denotes species, j denotes base point of trapezoid
c
      a11 = a1f(i,r0)
      x1  = 0.0
c
c     j loop is for the n radii in r. We save the integral at
c     each of these.  i points to the interval of ra that we are in.
c
      do 200 j=1,n
c
  110 if (r(j) .gt. ra(i+1))  go to 120
c
c     r(j) is within this interval of ra
c
      a12    = a1f(i,r(j))
      x2     = x(r(j))
      a1int  = a1int + (a11+a12)*(x2-x1)
      at1(j) = 0.5*a1int
      a11    = a12
      x1     = x2
      go to 200
c
c          r(j) is not in this interval, move to next i
c
  120 i     = i + 1
      a12   = a1(i)
      x2    = x(ra(i))
      a1int = a1int + (a11+a12)*(x2-x1)
      a11   = a12
      x1    = x2
      if (i .lt. nra)  go to 110
      if (j .lt. n)  call nuerr(4)
c
c     last r point was (rounded slightly) above ra(nra)
c
      at1(n) = 0.5 * a1int
  200 continue
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine atv2 (r, n, r0, at1, at2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     ATV2 calculates the value of the At integral from
c     r0 = rholt*SIN (theta) out to each of the nr radii specified
c     in r. This saves on redundant integration along the same
c     path. The integration is done by the trapezoidal rule.
c     Both species are done at once to save on overhead.
c
      include '../inc/neucg.m'
      dimension       r(*), at1(*), at2(*)
      real*8          singloid ! single precision variable
      common /nuatbl/ odelr, ra(nxmax), a1(nxmax), a2(nxmax), nra, nram1
c
c     These statement functions interpolate the A1 and A2 rates
c     to the desired point.  They assume an equally spaced grid with
c     ra(1) = 0 and ra(2) = delta r, as produced by the routine Aspace.
c     (odelr = 1/delta r)
c     Note that x is a leg of the right triangle with hypotenuse
c     rho and the other leg r0, thus following the neutral's path.
c
      a1f(i,rdum  ) = a1(i) + (a1(i+1)-a1(i)) * (rdum-ra(i))*odelr
      a2f(i,rdummy) = a2(i) + (a2(i+1)-a2(i)) * (rdummy-ra(i))*odelr
      x(dumr)       = SQRT (dumr*dumr - r0sq)
c
c     We start at r0
c
      r0sq     = r0 * r0
      a1int    = 0.0
      a2int    = 0.0
      singloid = r0 * odelr
      i        = MIN0 (IFIX (singloid) + 1, nram1)
c
c     aij: i denotes species, j denotes base point of trapezoid
c
      a11 = a1f(i,r0)
      a21 = a2f(i,r0)
      x1  = 0.0
c
c     j loop is for the n radii in r. We save the integral at
c     each of these.  i points to the interval of ra that we are in.
c
      do 200 j=1,n
  110 if (r(j) .gt. ra(i+1))  go to 120
c
c     r(j) is within this interval of ra
c
      a12 = a1f(i,r(j))
      a22 = a2f(i,r(j))
      x2 = x(r(j))
      dx = x2-x1
      a1int = a1int + (a11+a12)*dx
      a2int = a2int + (a21+a22)*dx
      at1(j) = 0.5*a1int
      at2(j) = 0.5*a2int
      a11 = a12
      a21 = a22
      x1 = X2
      go to 200
  120 continue
c
c          r(j) is not in this interval, move to next i
c
      i = i+1
      a12 = a1(i)
      a22 = a2(i)
      x2 = x(ra(i))
      dx = x2-x1
      a1int = a1int + (a11+a12)*dx
      a2int = a2int + (a21+a22)*dx
      a11 = a12
      a21 = a22
      x1 = X2
      if (i .lt. nra)  go to 110
      if (j .lt. n  )  call nuerr(4)
c
c     last r point was (rounded slightly) above ra(nra)
c
      at1(n) = 0.5*a1int
      at2(n) = 0.5*a2int
  200 continue
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine kgrid (iseg, rr, amr, rk, wk)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     Kgrid sets up the Krylov and Paltsev quadrature base points
c     and weights for each point rr = r(i). The weights wk include
c     a factor rho, but not the normalization for interval length.
c
      include '../inc/neucg.m'
c
      parameter      (nkint = 10)
      dimension  rk(*), wk(*)
      dimension  x(nkint), w(nkint)
      data       x/0.101394196e-1, 0.592766049e-1, 0.147205025   ,
     .             0.266893349   , 0.408097323   , 0.558341000   ,
     .             0.704006629   , 0.831548182   , 0.928781769   ,
     .             0.986191429  /,
     .           w/0.284355910e-1, 0.693868683e-1, 0.105280045   ,
     .             0.132373448   , 0.147936959   , 0.150260429   ,
     .             0.138791499   , 0.114239061   , 0.786206825e-1,
     .             0.352595508e-1/
c
      if (iseg .eq. 2)  go to 110
c
c     iseg = 1 for 0 to rr segment
c
      do i=1,nkint
        rk(i) = rr    * (1.0 - x(i))
        wk(i) = rk(i) * w(i)
      end do
      return
c
c     iseg = 2 for rr to rminor segment
c
  110 do i=1,nkint
        rk(i) = rr + amr * x(i)
        wk(i) = rk(i) * w(i)
      end do
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine kintth (plt, pgt, ir, ti1, ti1w, ti2, ti2w)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c KINTTH does the integral over theta for this pair of plt,pgt.
c The radius for vth is r(ir).  We use the fast approximations
c to the At functions, as they yield sufficient accuracy for
c smooth plasma profiles
c
      include '../inc/neucg.m'
      logical         twospc, onespc
      common /nuspcs/ twospc, onespc
      common /nuquad/ r(nr),wr(nr),rwr(nr),sint(nt),wt(nt)
      common /nuovth/ ovth1(nr),ovth2(nr)
c
      ti1  = 0.0
      ti1w = 0.0
      ti2  = 0.0
      ti2w = 0.0
      pgt2 = pgt * pgt
c
      do 100 k=1,nt
      r0  = plt * sint(k)
      wtx = wt(k) / SQRT (pgt2 - r0*r0)
      if (onespc)  call fat1(plt,pgt,r0,ap1,am1)
      if (twospc)  call fat2(plt,pgt,r0,ap1,am1,ap2,am2)
      apovr = ap1*ovth1(ir)
      amovr = am1*ovth1(ir)
      call g0g2(apovr,g0ap,g2ap)
      call g0g2(amovr,g0am,g2am)
      ti1   = ti1 + wtx*(g0ap + g0am)
      ti1w  = ti1w + wtx*(g2ap + g2am)
      if (onespc)  go to 100
      apovr = ap2*ovth2(ir)
      amovr = am2*ovth2(ir)
      call g0g2(apovr,g0ap,g2ap)
      call g0g2(amovr,g0am,g2am)
      ti2   = ti2 + wtx*(g0ap + g0am)
      ti2w  = ti2w + wtx*(g2ap + g2am)
  100 continue
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine fat1 (plt, pgt, plst, ap1, am1)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     FAT1 performs a fast approximate integration for the At
c     functions.  In computing the singular part of the kernel, with a
c     smooth plasma profile, this accuracy is sufficient. Note that we
c     actually return the functions A+ and A-. The short A- path often
c     dominates, and this way it is treated more accurately. The
c     longer A+ path is often strongly attenuated, and so it requires
c     a less accurate evaluation. The integration is done with one
c     Simpson's panel, unless both end points lie in the same interval
c     of ra, in which case we use a trapezoidal rule.
c
      include '../inc/neucg.m'
      common /nuatbl/ odelr, ra(nxmax), a1(nxmax), a2(nxmax), nra, nram1
c
c     This statement function interpolates the A1 rate to the desired
c     point. It assumes an equally spaced grid with ra(1) = 0 and
c     ra(2) = delta r, as produced by the routine Aspace. (odelr =
c     1/delta r) Note that x is a leg of the right triangle with
c     hypotenuse rho and the other leg r0, thus following the
c     neutral's path.
c
      a1f (i, r)  = a1(i) + (a1(i+1)-a1(i)) * (r-ra(i))*odelr
      x  (pdum  ) = SQRT (pdum*pdum - plst2)
      p  (xdum  ) = SQRT (xdum*xdum + plst2)
****  ia (dummyp) = MIN0 (IFIX (dummyp * odelr) + 1, nram1)
c
c     Here we use the convention plt = rholt, pgt = rhogt, and
c     plst = rholt * SIN (theta).  p = rho is the radius out to x.
c
      plst2  = plst * plst
      iplt   = ia (plt , odelr, nram1)
      ipgt   = ia (pgt , odelr, nram1)
      iplst  = ia (plst, odelr, nram1)
      a1plt  = a1f(iplt,plt)
      a1pgt  = a1f(ipgt,pgt)
      a1plst = a1f(iplst,plst)
      xplt   = x(plt)
      xpgt   = x(pgt)
c
c          note x(plst) = 0.0
c
c          Do the inner segment of the path (from plst to plt).
c          Note that A+ = 2 Ainner + A-
c
      if (iplst .ne. iplt)  go to 100
      ap1 = (a1plst + a1plt)*xplt
      go to 110
c
c          get midpoint for Simpson's panel
c
  100 xmid  = xplt*0.5
      pmid  = p(xmid)
      apmid = a1f (ia (pmid, odelr, nram1), pmid)
      ap1   = (a1plst + 4.0 * apmid + a1plt)*xplt*0.33333333
c
c          Now do the outer segment of the path (= A-) from plt to pgt
c
  110 if (iplt .ne. ipgt)  go to 120
      am1 = 0.5*(a1plt + a1pgt)*(xpgt - xplt)
      go to 200
c
c          get midpoint for Simpson's panel
c
  120 xmid  = (xplt + xpgt)*0.5
      pmid  = p(xmid)
      apmid = a1f (ia (pmid, odelr, nram1), pmid)
      am1   = (a1plt + 4.0 * apmid + a1pgt)*(xpgt - xplt)*0.16666667
c
c          complete A+
c
  200 ap1 = ap1 + am1
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine fat2 (plt, pgt, plst, ap1, am1, ap2, am2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     FAT2 performs a fast approximate integration for the At
c     functions.  In computing the singular part of the kernel,
c     with a smooth plasma profile, this accuracy is sufficient.
c     Note that we actually return the functions A+ and A-.
c     The short A- path often dominates, and this way it is treated
c     more accurately. The longer A+ path is often strongly
c     attenuated, and so it requires a less accurate evaluation.
c     The integration is done with one Simpson's panel, unless both
c     end points lie in the same interval of ra, in which case
c     we use a trapezoidal rule.
c     Both species are done at the same time to save on overhead.
c
      include '../inc/neucg.m'
c
      common /nuatbl/ odelr, ra(nxmax), a1(nxmax), a2(nxmax), nra, nram1
c
c     These statement functions interpolate the A1 and A2 rates
c     to the desired point.  They assume an equally spaced grid with
c     ra(1) = 0 and ra(2) = delta r, as produced by the routine Aspace.
c     (odelr =1/delta r)
c     Note that x is a leg of the right triangle with hypotenuse
c     rho and the other leg r0, thus following the neutral's path.
c
      a1f (i, r ) = a1(i) + (a1(i+1)-a1(i)) * (r-ra(i))*odelr
      a2f (i, r ) = a2(i) + (a2(i+1)-a2(i)) * (r-ra(i))*odelr
      x  (dummyp) = SQRT (dummyp*dummyp - plst2)
      p  (xdum  ) = SQRT (xdum  *xdum   + plst2)
****  ia (pdum  ) = MIN0 (IFIX (pdum * odelr) + 1, nram1)
c
c     Here we use the convention plt = rholt, pgt = rhogt, and
c     plst = rholt*SIN (theta).  p = rho is the radius out to x.
c
      plst2  = plst*plst
      iplt   = ia (plt , odelr, nram1)
      ipgt   = ia (pgt , odelr, nram1)
      iplst  = ia (plst, odelr, nram1)
      a1plt  = a1f(iplt,plt)
      a1pgt  = a1f(ipgt,pgt)
      a1plst = a1f(iplst,plst)
      a2plt  = a2f(iplt,plt)
      a2pgt  = a2f(ipgt,pgt)
      a2plst = a2f(iplst,plst)
      xplt   = x(plt)
      xpgt   = x(pgt)
c
c          note x(plst) = 0.0
c
c          Do the inner segment of the path (from plst to plt).
c          Note that A+ = 2 Ainner + A-
c
      if (iplst .eq. iplt) then
        ap1 = (a1plst + a1plt) * xplt
        ap2 = (a2plst + a2plt) * xplt
        go to 110
      end if
c
c     get midpoint for Simpson's panel
c
      xmid  = xplt * 0.5
      pmid  = p  (xmid)
      imid  = ia (pmid, odelr, nram1)
      apmid = a1f(imid,pmid)
      xw    = xplt * 0.33333333
      ap1   = (a1plst + 4.0 * apmid + a1plt) * xw
      apmid = a2f(imid,pmid)
      ap2   = (a2plst + 4.0 * apmid + a2plt) * xw
c
c     Now do the outer segment of the path (= A-) from plt to pgt
c
  110 if (iplt .ne. ipgt)  go to 120
      am1 = 0.5*(a1plt + a1pgt)*(xpgt - xplt)
      am2 = 0.5*(a2plt + a2pgt)*(xpgt - xplt)
      go to 200
c
c     get midpoint for Simpson's panel
c
  120 xmid  = (xplt + xpgt)*0.5
      pmid  = p  (xmid)
      imid  = ia (pmid, odelr, nram1)
      apmid = a1f(imid,pmid)
      xw    = (xpgt - xplt)*0.16666667
      am1   = (a1plt + 4.0 * apmid + a1pgt)*xw
      apmid = a2f(imid,pmid)
      am2   = (a2plt + 4.0 * apmid + a2pgt)*xw
c
c     complete A+
c
  200 ap1 = ap1 + am1
      ap2 = ap2 + am2
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      integer function ia (alpha, beta, iota)
c
      implicit none
c
c --- replacement for statement functions (below) in FAT1 and FAT2:
c --- ia (dummyp) = MIN0 (IFIX (dummyp * odelr) + 1, nram1) --- in FAT1
c --- ia (pdum  ) = MIN0 (IFIX (pdum   * odelr) + 1, nram1) --- in FAT2
c
      integer   iota
      real*8    singloid ! single precision variable
      real*8    alpha, beta
c
      singloid = alpha * beta
      ia       = MIN0 (IFIX (singloid) + 1, iota)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine icsicu1 (x, y, nx, bpar, c, ic, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- modified version of IMSL subroutine named ICSICU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - interpolatory approximation by cubic splines
c                           with arbitrary second derivative end
c                           conditions.
c
c   usage               - call icsicu1 (x, y, nx, bpar, c, ic, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx. (input) x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points.
c                           (input)
c                nx     - number of elements in x and y. (input) nx
c                           must be .ge. 2.
c                bpar   - vector of length 4 containing the end
c                           condition parameters. (input)
c                           2.0*spp(1)+bpar(1)*spp(2) = bpar(2),
c                           bpar(3)*spp(nx-1)+2.0*spp(nx) = bpar(4),
c                           where spp(i) = second derivative of the
c                           cubic spline function s evaluated at x(i).
c                c      - spline coefficients. (output) c is an nx-1 by
c                           3 matrix. the value of the spline
c                           approximation at t is
c                           s(t) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
c                           where x(i) .le. t .lt. x(i+1) and
c                           d = t-x(i).
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement in
c                           the calling program. (input)
c                ier    - error parameter. (output)
c                         terminal error
c                           ier = 129, ic is less than nx-1
c                           ier = 130, nx is less than 2.
c                           ier = 131, input abscissa are not ordered
c                             so that x(1) .lt. x(2) ... .lt. x(nx).
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,ier
      real*8             x(nx),y(nx),bpar(4),c(ic,3)
c
c     specifications for local variables
c
      integer            i,j,nxm1
      real*8             dx,dxj,dxjp1,dxp,dyj,dyjp1,half,one,pj,
     .                   six,sixi,two,yppa,yppb,zero
      equivalence        (dxj,yppb),(pj,sixi),(dxjp1,yppa)
      data               zero/0.0/,half/0.5/,one/1.0/,
     .                   two/2.0/,six/6.0/
c
      ier = 0
c
c     check error conditions
c
      nxm1 = nx-1
      if (ic .lt. nxm1)  go to 30
      if (nx .lt. 2   )  go to 35
      if (nx .eq. 2   )  go to 10
c
c     compute coefficients and right hand side of the tridiagonal
c     system defining the second derivatives of the spline interpolant for (x,y)
c
c     c(j,1) = lambda(j)
c     c(j,2) = mu(j)
c     c(j,3) = d(j)
c
      dxj = x(2)-x(1)
      if (dxj .le. zero)  go to 40
      dyj = y(2)-y(1)
      do 5 j=2,nxm1
         dxjp1 = x(j+1)-x(j)
         if (dxjp1 .le. zero)  go to 40
         dyjp1 = y(j+1)-y(j)
         dxp = dxj+dxjp1
         c(j,1) = dxjp1/dxp
         c(j,2) = one-c(j,1)
         c(j,3) = six*(dyjp1/dxjp1-dyj/dxj)/dxp
         dxj = dxjp1
         dyj = dyjp1
    5 continue
c
c     factor the tridiagonal matrix and solve for u
c
c     c(j,2)  = u(j)
c     c(j,1)  = q(j)
c     bpar(1) = lambda(1)
c     bpar(2) = d(1)
c     bpar(3) = mu(nx)
c     bpar(4) = d(nx)
c
   10 c(1,1) = -bpar(1)*half
      c(1,2) = bpar(2)*half
      if (nx .eq. 2)  go to 20
      do 15 j=2,nxm1
         pj = c(j,2)*c(j-1,1)+two
         c(j,1) = -c(j,1)/pj
         c(j,2) = (c(j,3)-c(j,2)*c(j-1,2))/pj
   15 continue
c
c     solve for cubic coefficients of spline interpolant
c     c(j,1), c(j,2), and c(j,3)
c
   20 yppb = (bpar(4)-bpar(3)*c(nxm1,2))/(bpar(3)*c(nxm1,1)+two)
      sixi = one/six
      do 25 i=1,nxm1
         j = nx-i
         yppa = c(j,1)*yppb+c(j,2)
         dx = x(j+1)-x(j)
         c(j,3) = sixi*(yppb-yppa)/dx
         c(j,2) = half*yppa
         c(j,1) = (y(j+1)-y(j))/dx-(c(j,2)+c(j,3)*dx)*dx
         yppb = yppa
   25 continue
      go to 9005
   30 ier = 129
      go to 9000
   35 ier = 130
      go to 9000
   40 ier = 131
c
 9000 continue
c     call uertst1 (ier, 'icsicu1')
 9005 return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine icsevu1 (x, y, nx, c, ic, u, s, m, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- modified version of IMSL subroutine named ICSEVU
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - evaluation of a cubic spline
c
c   usage               - call icsevu1 (x, y, nx, c, ic, u, s, m, ier)
c
c   arguments    x      - vector of length nx containing the abscissae
c                           of the nx data points (x(i),y(i)) i = 1,...,
c                           nx (input). x must be ordered so that
c                           x(i) .lt. x(i+1).
c                y      - vector of length nx containing the ordinates
c                           (or function values) of the nx data points
c                           (input).
c                nx     - number of elements in x and y (input).
c                           nx must be .ge. 2.
c                c      - spline coefficients (input). c is an nx-1 by
c                           3 matrix.
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement
c                           in the calling program (input).
c                           ic must be .ge. nx-1
c                u      - vector of length m containing the abscissae
c                           of the m points at which the cubic spline
c                           is to be evaluated (input).
c                s      - vector of length m (output).
c                           the value of the spline approximation at
c                           u(i) is
c                           s(i) = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
c                           where x(j) .le. u(i) .lt. x(j+1) and
c                           d = u(i)-x(j).
c                m      - number of elements in u and s (input).
c                ier    - error parameter (output).
c                         warning error
c                           ier = 33, u(i) is less than x(1).
c                           ier = 34, u(i) is greater than x(nx).
c
c                           ********************************************
c                           output of warning errors has been suppressed
c                           ********************************************
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   remarks  1.  the routine assumes that the abscissae of the nx
c                data points are ordered such that x(i) is less than
c                x(i+1) for i = 1,...,nx-1. no check of this condition
c                is made in the routine. unordered abscissae will cause
c                the algorithm to produce incorrect results.
c            2.  the routine generates two warning errors. one error
c                occurs if u(i) is less than x(1), for some i in the
c                the interval (1,m) inclusively. the other error occurs
c                if u(i) is greater than x(nx), for some i in the
c                interval (1,m) inclusively.
c            3.  the ordinate y(nx) is not used by the routine. for
c                u(k) .gt. x(nx-1), the value of the spline, s(k), is
c                given by
c                 s(k) = ((c(nx-1,3)*d+c(nx-1,2))*d+c(nx-1,1))*d+y(nx-1)
c                where d = u(k)-x(nx-1).
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
c     specifications for arguments
c
      integer            nx,ic,m,ier
      real*8             x(nx),y(nx),c(ic,3),u(m),s(m)
c
c     specifications for local variables
c
      integer            i,jer,ker,nxm1,k
      real*8             d,dd,zero
      data               i/1/, zero/0.0/
c
c     first executable statement
c
      jer = 0
      ker = 0
      if (m .le. 0)  go to 9005
      nxm1 = nx-1
      if (i .gt. nxm1)  i = 1
c
c     evaluate spline at m points
c
      do 40 k=1,m
c
c        find the proper interval
c
         d = u(k)-x(i)
         if (d) 5, 25, 15
    5    if (i .eq. 1)  go to 30
         i = i-1
         d = u(k)-x(i)
         if (d) 5, 25, 20
   10    i = i+1
         d = dd
   15    if (i .ge. nx)  go to 35
         dd = u(k)-x(i+1)
         if (dd .ge. zero)  go to 10
         if ( d .eq. zero)  go to 25
c
c        perform evaluation
c
   20    s(k) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+y(i)
         go to 40
   25    s(k) = y(i)
         go to 40
c
c        u(k) < x(1)
c
   30    jer = 33
         go to 20
c
c        u(k) > x(nx)
c
   35    if (dd .gt. zero)  ker = 34
         d = u(k) - x(nxm1)
         i = nxm1
         go to 20
c
   40 continue
c
      ier = MAX0 (jer, ker)
c
c***  if (jer .gt. 0)  call uertst1 (jer, 'icsevu1')
c***  if (ker .gt. 0)  call uertst1 (ker, 'icsevu1')
c
 9005 return
c
      end
