
      subroutine pelabl (rps,tes,dens,xrhor,rdot1,rdot2,difrdt,newcel)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- pelabl calculates the rates of recession of the pellet surface
c --- from two equations, using the cloud thickness as the independent parameter
c --- references:
c --- s.l.milora,ornl/tm-8616 (1983)
c --- last revision: 7/83 w.a.houlberg and s.e.attenberger ornl
c --- calculated parameters:
c --- rdot1-dr/dt from energy balance at pellet surface-(cm/s)
c --- rdot2-dr/dt from balance in cloud-(cm/s)
c --- difrdt-ratio of difference to average of the two dr/dt solutions
c --- input parameters:
c --- xrhor-cloud thickness / rhosrp-(dimensionless)
c --- other parameters:
c --- xmp-molecular mass of pellet species-(g/molecule)
c --- hvap-pellet heat of evaporation-(eV/g)
c --- other comments:
c ----------------------------------------------------------------------
c
      common /compll/ rpel, amup, denm, rhosrp, xrhoro, qeo, fqes, fqe,
     .                eion, ebg(18), denbg(17), ebgo(17), qbgo(17), ab,
     .                vcb, qbo, fqb, eag(21), denag(20), eago(20),
     .                qago(20), aa, vca, qao, fqa,
     .                nag, nbg
c
      data            gam/1.4/, evap/0.01/
c
      pi   = ACOS (-1.0)
      xmp  = 2.0 * amup*1.6726e-24
      hvap = evap / xmp
c
c --- set pellet surface area exposed to each flux.
c
      areae = 2.0 * pi * rps**2
      areaf = 2.0 * areae
      areat = areaf
c
c --- get the electron and ion heat fluxes to cloud and pellet surface.
c
      call pelqe (tes,dens,xmp,rhosrp,xrhor,qeo,fqes,fqe)
      denbgt = ssumpl(nbg,denbg,1)
      if (denbgt .le. 0.0)  go to 10
      call pelqf (xmp,ebg,denbg,nbg,ab,vcb,rhosrp,xrhor,ebgo,qbgo,
     .            qbo,fqb,1,newcel)
   10 denagt = ssumpl(nag,denag,1)
      if (denagt .le. 0.0)  go to 20
      call pelqf (xmp,eag,denag,nag,aa,vca,rhosrp,xrhor,eago,qago,
     .            qao,fqa,2,newcel)
c
c --- calculate the pellet surface erosion rate
c
   20 qtp    =  qeo*fqe*areae+(qbo*fqb+qao*fqa)*areaf
      rdot1  = -qtp*rpel/(areat*hvap*rhosrp)
      qtc    =  qeo*(1.0-fqe)+qbo*(1.0-fqb)+qao*(1.0-fqa)
      q      =  qtc/(rhosrp*xrhor)
      xr     =  rps/rpel
      rdot2  = -1.25*xrhor/xr*(1.602e-12*q*rps*(gam-1.0)/2.0)**(1.0/3.0)
      difrdt =  2.0*(rdot2-rdot1)/(rdot2+rdot1)
      return
c
      end

      subroutine pelfil (ndat, nin, versid)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- create pellet plot data file, and copy inone file text to it
c
      character*8  title(10)
      character *(*) versid
c
      call DESTROY ('peldat')
      open   (unit = ndat, file = 'peldat', status = 'UNKNOWN')
      write  (ndat, 8010) versid
 8010 format (' ******** ', a, ' ******** (33x65)' ///)
      rewind (unit = nin)
c
 2000 read   (nin , '(    9a8)', end=2020)  (title(i), i=1,9)
      write  (ndat, '(1x, 9a8)'          )  (title(i), i=1,9)
      go to 2000
c
 2020 write  (ndat, '(''stop'')')
      rewind (unit = nin)
      return
c
      end

      subroutine pellet (kprint,nprt,ndat,amupx,rpelx,vpel,n,kj,dvol,nl,
     .                   mapil,rchord,den,te,abx,nbgx,ebgx,vcbx,denbig,
     .                   aax,nagx,eagx,vcax,denaig,pden)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- PELLET calculates the particle deposition profile for frozen
c --- hydrogenic pellets injected at any angle in the plasma midplane
c --- of an axisymmetric toroidal system.
c --- references:
c --- s.l.milora,ornl/tm-8616 (1983).
c --- w.a.houlberg,m.a.iskra,h.c.howe,s.e.attenberger,ornl/tm-6549
c --- (1979).
c --- last revision: 7/83 w.a.houlberg and s.e.attenberger ornl.
c --- calculated parameters:
c --- pden(i)-increase in plasma density in region i-(/cm**3).
c --- input parameters:
c --- kprint-option for detail of printout.
c ---       = 0 errors in input and solution.
c ---       = 1 above plus input values and options.
c ---       = 2 above plus radial input and calculated parameters.
c ---       = 3 above plus parameters along pellet path.
c --- amupx-atomic mass of pellet atoms-(1<amup<3).
c --- rpelx-initial pellet radius-(cm).
c --- vpel-pellet velocity-(cm/sec).
c --- n-number of radial regions in plasma.
c --- dvol(i)-volume of plasma region i-(cm**3).
c --- nl-number of segments along pellet path.
c --- mapil(l)-plasma region i for pellet segment l.
c --- rchord(l)-path distance to start of segment l-(cm).
c ---           (rchord(nl+1) is the far wall and should be non-zero.)
c --- den(i)-electron density in plasma zone i-(/cm**3).
c --- te(i)-electron temperature in plasma zone i-(ev).
c --- abx-atomic mass number of fast h ions-(dimensionless).
c --- nbgx-number of fast h ion energy groups.
c --- ebgx(jg)-fast h ion energy at group boundary jg-(ev).
c ---         -ebgx(jg) > ebgx(jg+1).
c --- vcbx(i)-critical velocity in region i for fast h ions -(cm/s).
c --- denbig(i,jg)-fast h density, interval jg and region i-(/cm**3).
c --- aax-atomic mass number of fast he ions-(dimensionless).
c --- nagx-number of fast he ion energy groups.
c --- eagx(jg)-fast he ion energy at group boundary jg-(ev).
c ---         -eagx(jg) > eagx(jg+1).
c --- vcax(i)-critical velocity in region i for fast he ions -(cm/s).
c --- denaig(i,jg)-fast he density, interval jg and region i-(/cm**3).
c --- other parameters:
c --- denm-molecular density of hydrogenic ice-(/cm**3).
c ----------------------------------------------------------------------
c
      common /compll/ rpel, amup, denm, rhosrp, xrhoro, qeo, fqes, fqe,
     .                eion, ebg(18), denbg(17), ebgo(17), qbgo(17), ab,
     .                vcb, qbo, fqb, eag(21), denag(20), eago(20),
     .                qago(20), aa, vca, qao, fqa,
     .                nag, nbg
c
      common /units/  nout
      dimension       dvol(*),den(*),te(*),vcbx(*),vcax(*),pden(*)
      dimension       ebgx(18),denbig(kj,17)
      dimension       eagx(2),denaig(kj,*)
      dimension       rchord(*),mapil(*)
c
      pi     = ACOS (-1.0)
      eion   = 32.6
      nout   = 23
      amup   = amupx
      denm   = -8.6857e20*amup**2+6.3023e21*amup+2.1200e22
      rpel   = rpelx
      rhosrp = (2.0*amup*denm*1.6726e-24)*rpel
      ab     = abx
      nbg    = nbgx
      do jg=1,nbg+1
        ebg(jg) = ebgx(jg)
      end do
      aa  = aax
      nag = nagx
      do jg=1,nag+1
        eag(jg) = eagx(jg)
      end do
      do i=1,n
        pden(i) = 0.0
      end do
      nbgp1 = nbg+1
      t     = 0.0
      rp    = rpel
      xrhoro = 0.0
      qbo    = 0.0
      fqb    = 0.0
      qao    = 0.0
      fqa    = 0.0
      if (kprint .lt. 1)  go to 40
c
c --- optional output
c
      write (nprt, 1010)
      write (nprt, 1020) kprint,n,amup,rpel,vpel
      if (kprint .lt. 2)  go to 40
      write (nprt, 1030)
c
c --- follow the pellet path and determine the ablation rate
c
   40 do l=1,nl
        if (rp .le. 0.0)  go to 80
        i      = mapil(l)
        denold = den(i)+pden(i)
        dennew = denold
        teold  = (den(i)*te(i)-pden(i)*(2.0/3.0)*eion)/denold
        tenew  = teold
        vcb    = vcbx(i)
        do jg=1,nbg
          denbg(jg) = denbig(i,jg)
        end do
        vca = vcax(i)
        do jg=1,nag
          denag(jg) = denaig(i,jg)
        end do
        dt = (rchord(l+1)-rchord(l))/vpel
        call pelrk4 (dt, t, rp, tenew, dennew, dvol(i))
        if (rp .lt. 0.0)  rp = 0.0
        srcp    = dennew-denold
        src     = srcp*dvol(i)
        pden(i) = pden(i)+srcp
        if (kprint .gt. 2)
     .    write (nprt, 1040) l,mapil(l),rchord(l),teold,denold,rp,
     .                       src,srcp,qeo,fqe,fqes,qbo,fqb,xrhoro
        write  (ndat, 1070)  rchord(l),src,srcp,qeo,fqe,qbo,fqb
 1070   format (7(E16.8))
      end do
c
   80 dps = -1.0
      write  (ndat, 1071) dps
 1071 format (e16.8)
      write  (ndat, 1072) amup,rpel,vpel,n
 1072 format (3(e16.8),i3)
c
      if (kprint .lt. 2)  return
      write (nprt, 1050)
c
      do i=1,n
        dsdn = 0.0
        if (den(i) .gt. 0.0)  dsdn = pden(i)/den(i)
        write  (nprt, 1060)  i, te(i), den(i), dvol(i), pden(i), dsdn
        write  (ndat, 1073)  te(i), den(i), pden(i)
 1073   format (3(e16.8))
      end do
c
 1010 format (1h1 / 30x, '*** summary from subroutine PELLET ***' //)
 1020 format (20x, 'kprint -- print option               = ',   i5   /
     .        20x, 'n      -- number of mesh points      = ',   i5   /
     .        20x, 'amup   -- ave atomic mass of pellet  = ',  f10.2 /
     .        20x, 'rpel   -- initial pellet radius (cm) = ',  f10.4 /
     .        20x, 'vpel   -- pellet velocity (cm/sec)   =', 1pe10.3 )
 1030 format (1h1 / 30x, '*** summary along pellet path ***' //
     . ' path rad  path dist elec temp  density   pel rad   source',
     . '    source     qeo        fqe       fqes      qbo       fqb',
     . '     rhorhat' /
     . ' cell cell   (cm)      (eV)     (cm-3)      (cm)',
     . '              (cm-3)   (eV*cm-2)                     (eV*cm-2)')
 1040 format (2(2x,i3),2(1x,0pf9.2),10(1x,1pe9.2))
 1050 format (1h1 / 25x,41h*** summary of radial pellet profiles *** //
     .3x,'radial   elec temp  plasma den    volume      source     so',
     . 'urce',//
     .3x,'node      (ev)       (cm-3)      (cm+3)      (cm-3)     de',
     . 'nsity',/)
 1060 format (4x,i3,4x,4(1pe10.3,2x),0pf10.4)
      return
c
      end

      subroutine pelqe (te, den, xmp, rhosrp, xrhor, qeo, fqes, fqe)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- pelqe calculates the electron heat flux incident on the pellet
c --- cloud and the heat flux attenuation factor in the cloud.
c --- references:
c --- s.l.milora,ornl/tm-8616 (1983).
c --- last revision: 2/84 w.a.houlberg and s.e.attenberger ornl.
c --- calculated parameters:
c --- qeo-electron heat flux incident on pellet cloud-(ev/cm**2).
c --- fqes = qes/qeo-heat flux attenuation factor-(dimensionless).
c --- fqe = qep/qeo-heat flux attenuation factor-(dimensionless).
c --- input parameters:
c --- te-electron temperature in plasma-(ev).
c --- den-electron density in plasma-(/cm**3).
c --- xmp-molecular mass of pellet species-(g/molecule).
c --- rhosrp-pellet radius times solid pellet density-(g/cm**2).
c --- xrhor-cloud thickness/rhosrp-(dimensionless).
c --- other parameters:
c --- es-energy below which elastic scattering is included-(ev).
c --- qes-electron heat flux at energy es = 20 ev-(ev/cm**2).
c --- qep-electron heat flux at pellet surface-(ev/cm**2).
c --- alfe-cross-section for elastic scattering-(cm**2).
c --- ce-mean electron thermal speed in plasma-(cm/s).
c --- cjeo-random electron particle flux in plasma-(/cm**2/s).
c --- neg-number of electron groups
c ----------------------------------------------------------------------
c
      dimension we(1)
      data      a1e/2.35e14/, a2e/4.0e11/
      data      c/2.0/, es/100.0/, alfe/1.8e-16/
****  dimension we(5)
****  data we/3.58220,2.39789,1.76644,1.25078,0.50209/, neg/5/
      data      we/1.0/, neg/1/
c
      pi   = ACOS (-1.0)
      wq   = 1.0 / neg
      c1   = a1e/a2e
      c2   = 2.0 * rhosrp*xrhor/(xmp*a2e)
      eo   = 1.5*te
      ce   = SQRT ((8.0/pi)*te*1.602e-12/9.105e-28)
      cjeo = den*ce/4.0
      qeo  = (4.0/3.0)*cjeo*eo
      qeog = qeo*wq
      qes  = (4.0/3.0)*cjeo*es
      qep  = 0.0
c
c --- set xi1 for all energies
c
      xi1 = rhosrp*xrhor/xmp
c
      do jg=1,neg
        eog  = eo*we(jg)
        qesg = qeog*es/eog
c
c --- set xi2 and incident flux to scattering regime for eog < es
c
        xi2 = 0.0
        qex = qeog
        if (eog .lt. es)  go to 10
c
c --- set xi2 and incident flux to scattering regime for eog > es
c
        xi2 = a1e*(eog-es)+a2e*(eog**2-es**2)/2.0
        qex = qesg
   10   xi  = alfe*(xi1-xi2)
c
c --- check whether epg > es
c
        if (xi .lt. 0.0)  go to 20
c
c --- epg < es
c
        if (xi .gt. 85.0) xi = 85.0
        qepg = qex*(c+1)/(c+EXP (xi))
        go to 30
c
c --- epg > es
c
   20   epg  = -c1 + SQRT ((c1+eog)**2-c2)
        qepg = qeog*epg/eog
   30   qep  = qep + qepg
      end do
c
      fqes = qep/qes
      fqe  = qep/qeo
      return
c
      end

      subroutine pelqf (xmp,efg,denfg,nfg,af,vcf,rhosrp,
     .                  xrhor,efgo,qfgo,qfo,fqf,izf,newcel)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- pelqf calculates the fast ion heat flux incident on the pellet
c --- cloud and the heat flux attenuation factor in the cloud for
c --- either hydrogenic or helium fast ions.
c --- references:
c --- s.l.milora,ornl/tm-8616 (1983).
c --- last revision: 7/83 w.a.houlberg and s.e.attenberger ornl.
c --- calculated parameters:
c --- fqf = qfp/qfo-heat flux attenuation factor-(dimensionless).
c --- calculated parameters (newcel = 1):
c --- efgo(jg)-average   energy in energy interval jg-(ev).
c --- qfgo(jg)-heat flux        in energy interval jg-(/cm**2/s).
c --- qfo-fast ion heat flux incident on pellet cloud-(ev/cm**2).
c --- input parameters:
c --- xmp-molecular mass of pellet species-(g/molecule).
c --- efg(jg)-fast ion energy at group boundary jg-(ev).
c ---        -efg(jg) > efg(jg+1).
c --- denfg(jg)-fast ion density in energy interval jg-(/cm**3).
c ---          -efg(jg) > e > efg(jg+1).
c --- nfg-number of fast energy groups.
c --- af-atomic mass number of fast ions-(dimensionless).
c --- vcf-critical velocity for classical thermalization-(cm/s).
c --- rhosrp-pellet radius times solid pellet density-(g/cm**2).
c --- xrhor-cloud thickness / rhosrp-(dimensionless).
c --- newcel-flag for reevaluating fluxes incident on cloud.
c ---       =0 use old values of qfo, efgo and qfgo from earlier call.
c ---       =1 reset qfo, efgo and qfgo for new cell.
c --- izf-fast ion charge-(dimensionless).
c ---   =1 fast hydrogen ions.
c ---   =2 fast helium ions.
c --- input parameters (newcel = 0):
c --- efgo(jg)-see above.
c --- qfgo(jg)-see above.
c --- qfo-see above.
c --- other parameters:
c --- a0h,a1h,a2h-parameters for fit to fast h+ energy loss in h2 gas.
c --- a1he,a2he-parameters for fit to fast he+ energy loss in h2 gas.
c --- qfp-fast ion heat flux at pellet surface-(ev/cm**2).
c ----------------------------------------------------------------------
c
      dimension efg(*), denfg(*), efgo(*), qfgo(*)
      data      a0h/1.0e-13/, a1h/978.2/, a2h/3.239/
      data      a1he/2.788e-16/, a2he/0.2968/
c
      if ((izf .lt. 1) .or. (izf .gt. 2))  return
      if (newcel .eq. 0)  go to 20
c
c --- calculate parameters incident on pellet cloud
c
      cjo = 0.0
      qfo = 0.0
      sq3 = SQRT (3.0)
c
c --- set velocities and integral expressions at low end of group.
c
      vl   = SQRT (2.0*efg(1)*1.602e-12/(1.6726e-24*af))
      xvl  = vl/vcf
      xvl2 = xvl**2
      xvl3 = xvl*xvl2
      tl1  =  LOG (1.0+xvl3)
      tl2  =  LOG ((1.0-xvl+xvl2)/(1.0 + xvl)**2)
      tl3  = ATAN ((2.0*xvl-1.0)/sq3)
c
c --- sum over energy groups
c
      do jg=1,nfg
c
c --- shift expressions to boundaries of next lower energy group.
c
        xvh  = xvl
        vl   = SQRT (2.0*efg(jg+1)*1.602e-12/(1.6726e-24*af))
        xvl  = vl/vcf
        xvh2 = xvl2
        xvl2 = xvl**2
        xvh3 = xvl3
        xvl3 = xvl*xvl2
        th1  = tl1
        tl1  = LOG (1.0+xvl3)
        th2  = tl2
        tl2  = LOG ((1.0-xvl+xvl2)/(1.0 + xvl)**2)
        th3  = tl3
        tl3  = ATAN ((2.0*xvl-1.0)/sq3)
c
c --- calculate mean group energy in plasma
c
        cong     = 3.0 / LOG ((1.0 + xvh3) / (1.0 + xvl3))
        efgo(jg) = 1.6726e-24*af*0.5*vcf**2*cong/1.602e-12
     .           * ((xvh2*0.5-(0.5*th2+sq3*th3)/3.0)
     .           - (xvl2*0.5-(0.5*tl2+sq3*tl3)/3.0))
c
c --- calculate mean group heat flux in plasma
c
        qfgo(jg) = denfg(jg)*1.6726e-24*af*vcf**3/(1.602e-12*24.0)
     .           * cong *((xvh3-th1)-(xvl3-tl1))
        qfo      = qfo+qfgo(jg)
      end do
c
c --- calculate parameters incident on the pellet
c
   20 qfp = 0.0
c
c --- check fast ion species
c
      if (izf .eq. 2)  go to 40
c
c --- fast hydrogenic ions
c
      b  = a1h * SQRT (af) / a2h
      bb = b**2
      arhmg = a0h*rhosrp*xrhor/xmp/a2h
      do jg=1,nfg
        egp = 0.0
        c   = 2.0 * SQRT (efgo(jg))/b+(efgo(jg)-arhmg)/bb
        if (c .gt. 0.0) egp = bb*(-1.0 + SQRT (1.0+c))**2
        qfp = qfp+qfgo(jg)*egp/efgo(jg)
      end do
      fqf = qfp/qfo
      return
c
c --- fast helium ions
c
   40 c1 = 1.0-a2he
      c2 = (4.0/af)**a2he*a1he*c1*rhosrp*xrhor/xmp
      do jg=1,nfg
        egp = 0.0
        efgot = efgo(jg)**c1
        if (efgot .gt. c2) egp = (efgot-c2)**(1.0/c1)
        qfp = qfp+qfgo(jg)*egp/efgo(jg)
      end do
      fqf = qfp/qfo
      return
c
      end

      subroutine pelrat (rps, tes, dens, rdot, newcel)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- pelrat determines the pellet ablation rate through iterating on
c --- the cloud thickness.
c --- references:
c --- s.l.milora,ornl/tm-8616 (1983).
c --- forsythe,malcolm,moler p.161.
c --- last revision: 7/83 w.a.houlberg and s.e.attenberger ornl.
c --- calculated parameters:
c --- rdot-change in pellet radius in time-(cm/sec).
c --- input parameters:
c --- newcel-option to recalculate fast ion fluxes in plasma.
c ---     =0 same plasma cell as previous call.
c ---     =1 new plasma cell.
c --- other comments:
c --- the zeroin procedure from forsythe,et al., is used to find the
c --- simultaneous solution to two equations for dr/dt which are both
c --- nonlinear functions of the cloud thickness and plasma parameters.
c ----------------------------------------------------------------------
c
      common /compll/ rpel, amup, denm, rhosrp, xrhoro, qeo, fqes, fqe,
     .                eion, ebg(18), denbg(17), ebgo(17), qbgo(17), ab,
     .                vcb, qbo, fqb, eag(21), denag(20), eago(20),
     .                qago(20), aa, vca, qao, fqa,
     .                nag, nbg
c
      data tolpr/0.02/
      data amin/1.0e-10/, bmax/1.0e-1/, xrange/3.0/
c
      one = 1.0
c
c --- check pellet size
c
      if ((rps .lt. (tolpr*rpel)) .or. (tes .lt. 1.0))  go to 130
c
c --- compute the relative machine precision
c
      eps = 1.0
   10 eps = eps/2.0
      test1 = 1.0 + eps
      if (test1 .gt. 1.0)  go to 10
c
c --- begin zeroin procedure
c
      if (xrhoro .gt. amin)  go to 30
c
c --- use entire range
c
   20 a = amin
      b = bmax
      call pelabl(rps,tes,dens,a,rdot1a,rdot2a,difa,newcel)
      newcel = 0
      call pelabl(rps,tes,dens,b,rdot1b,rdot2b,difb,newcel)
      if ((SIGN (one, difa) + SIGN (one, difb)) .ne. 0.0)  go to 1000
      go to 50
c
c --- evaluate at previous solution
c
   30 b = xrhoro
      call pelabl(rps,tes,dens,b,rdot1b,rdot2b,difb,newcel)
      newcel = 0
      if (ABS (difb) .lt. tolpr)  go to 120
      if (ABS (difb) .ge. (2.0-tolpr))  go to 20
      a = b
      difa = difb
      rdot2a = rdot2b
c
c --- evaluate at previous solution*xrange.
c
      b = xrhoro*xrange
      if (b .gt. bmax) b = bmax
      call pelabl(rps,tes,dens,b,rdot1b,rdot2b,difb,newcel)
      if ((SIGN (one, difa) + SIGN (one, difb)) .eq. 0.0)  go to 50
c
c --- check which side to try next
c
      if (ABS (difb) .gt. ABS (difa))  go to 40
c
c --- evaluate at bmax.
c
      a = b
      difa = difb
      rdot2a = rdot2b
      b = bmax
      call pelabl(rps,tes,dens,b,rdot1b,rdot2b,difb,newcel)
      if ((SIGN (one, difa) + SIGN (one, difb)) .ne. 0.0)  go to 1000
      go to 50
c
c --- evaluate at previous solution/xrange.
c
   40 b = a
      difb = difa
      rdot2b = rdot2a
      a = xrhoro/xrange
      if (a .lt. amin) a = amin
      call pelabl(rps,tes,dens,a,rdot1a,rdot2a,difa,newcel)
      if ((SIGN (one, difa) + SIGN (one, difb)) .eq. 0)  go to 50
c
c --- evaluate at amin.
c
      b = a
      difb = difa
      rdot2b = rdot2a
      a = amin
      call pelabl(rps,tes,dens,a,rdot1a,rdot2a,difa,newcel)
      if ((SIGN (one, difa) + SIGN (one, difb)) .ne. 0)  go to 1000
c
c --- begin step.
c
   50 c = a
      difc = difa
      rdot2c = rdot2a
      d = b-a
      e = d
   60 if (ABS (difc) .ge. ABS (difb))  go to 70
      a = b
      b = c
      c = a
      difa = difb
      rdot2a = rdot2b
      difb = difc
      rdot2b = rdot2c
      difc = difa
      rdot2c = rdot2a
c
c --- convergence test.
c
   70 tltst = 2.0 * eps * ABS (b)
      fm = 0.5*(c-b)
      if (ABS (fm) .le. tltst) write (6,1020) fm,tltst
      if (ABS (fm) .le. tltst)  go to 120
      if (ABS (difb) .lt. tolpr)  go to 120
c
c --- check if bisection is necessary.
c
      if (ABS (e) .lt. tltst)  go to 100
      if (ABS (difa) .le. ABS (difb))  go to 100
c
c --- check if quadratic interpolation is possible.
c
      if (a .ne. c)  go to 80
c
c --- linear interpolation.
c
      s = difb/difa
      p = 2.0 * fm*s
      q = 1.0 - s
      go to 90
c
c --- inverse quadratic interpolation.
c
   80 q = difa/difc
      r = difb/difc
      s = difb/difa
      p = s*(2.0*fm*q*(q-r)-(b-a)*(r-1.0))
      q = (q-1.0)*(r-1.0)*(s-1.0)
c
c --- adjust signs.
c
   90 if (p .gt. 0.0) q = -q
      p = ABS (p)
c
c --- check if interpolation is acceptable.
c
      if ((2.0*p) .ge. (3.0*fm*q - ABS (tltst*q)))  go to 100
      if (p .ge. ABS (0.5*e*q))  go to 100
      e = d
      d = p/q
      go to 110
c
c --- bisection.
c
  100 d = fm
      e = d
c
c --- complete step.
c
  110 a = b
      difa = difb
      rdot2a = rdot2b
      if (ABS (d) .gt. tltst) b = b+d
      if (ABS (d) .le. tltst) b = b + SIGN (tltst,fm)
      call pelabl(rps,tes,dens,b,rdot1b,rdot2b,difb,newcel)
      if ((difb*(difc / ABS (difc))) .gt. 0.0)  go to 50
      go to 60
c
c --- end of zeroin procedure
c
  120 xrhor = b
      xrhoro = xrhor
      rdot = rdot2b
      return
c
c --- default calculation - no ablation
c
  130 rdot = 0.0
      return
 1000 write  (6, 1010)  difa, difb
 1010 format (' solution out of range in PELRAT, difa,difb = ', 2e12.3)
 1020 format (' no convergence in PELRAT, fm, tltst = ',        2e12.3)
      call STOP ('subroutine PELRAT: lack of convergence', 62)
c
      end

      subroutine pelrk4 (dt, t, rp, te, den, dvol)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- pelrk4 is a fourth order Runge-Kutta integration routine which
c --- updates the pellet radius and plasma parameters for time dt.
c --- references:
c --- s.l.milora,ornl/tm-8616 (1983).
c --- last revision: 7/83 w.a.houlberg and s.e.attenberger ornl.
c --- calculated parameters (updated):
c --- t-time since injection-(s).
c --- rp-pellet radius-(cm).
c --- te-electron temperature-(ev).
c --- den-electron density-(/cm**3).
c --- input parameters:
c --- dt-time the pellet spends in this plasma region-(sec).
c --- dvol-volume of this plasma region-(cm**3).
c --- other comments:
c ----------------------------------------------------------------------
c
      common /compll/ rpel, amup, denm, rhosrp, xrhoro, qeo, fqes, fqe,
     .                eion, ebg(18), denbg(17), ebgo(17), qbgo(17), ab,
     .                vcb, qbo, fqb, eag(21), denag(20), eago(20),
     .                qago(20), aa, vca, qao, fqa,
     .                nag, nbg
c
      dimension dr(4)
      data      tolpte/0.1/, tolpr/0.02/
c
      s(r1,r2) = (4.0/3.0) * pi * (r2**3-r1**3)*2.0*denm/dvol
c
c --- initialize for one step through cell.
c
      pi     = ACOS (-1.0)
      newcel = 1
      tnew   = t+dt
      dts    = dt
      nstep  = 0
c
c --- initialize intermediate Runge-Kutta parameters.
c
   10 rps   = rp
      tes   = te
      dens  = den
      nstep = nstep+1
c
c --- advance pellet through time interval dts.
c
      do i=1,4
        call pelrat(rps,tes,dens,rdot,newcel)
        dr(i) = dts*rdot
        rps   = rp+dr(i)/2.0
        if (i .gt. 2)  rps = rp+dr(i)
        dden  = s(rps,rp)
        dnte  = dden*(2.0/3.0)*eion
c
c --- check perturbation on background plasma.
c
        if (i .gt. 1)  go to 20
        fte  = 2.0 * (dden+dnte/te)/(den+2.0*dden)
        if (fte .gt. 1.0) fte = 1.0
        fdts = 1.0
        if ( fte .ne. 0.0)  fdts = tolpte/fte
        if (fdts .lt. 1.0)  go to 40
   20   dens = den+dden
        tes  = (den*te-dnte)/dens
      end do
c
c --- completed step
c
      rps  = rp+(dr(1)+2.0*(dr(2)+dr(3))+dr(4))/6.0
      if (rps .lt. (tolpr*rpel))  rps = 0.0
      dden = s(rps,rp)
      dnte = dden*(2.0/3.0)*eion
      dens = den+dden
      tes  = (den*te-dnte)/dens
      den  = dens
      te   = tes
      rp   = rps
      t    = t+dts
      if ((t .ge. tnew) .or. (rp .eq. 0.0))  return
      if ((t+dts) .gt. tnew)  dts = tnew-t
      go to 10
c
c --- reduce time step size and try again.
c
   40 dts = 0.9*fdts*dts
      go to 10
c
      end
      subroutine propel

c ----------------------------------------------------------------------
c     propel         mar. 11, 1984       Dave S., Bob S.
c
c     ONETWO interface for Oak Ridge pellet fueling model.  the pellet
c     model needs electron density and temperature profiles, and fast
c     ion density profiles resolved into energy bins.  the model returns
c     the increase in electron density in each plasma shell due to the
c     complete ablation of a pellet.
c
c     other ONETWO inputs via common /pelcom/:
c     ipel   - pellet ion species index
c     pelrad - initial pellet radius (cm)
c     vpel   - pellet velocity (cm/sec)
c     nbgpel - number of fast beam ion energy bins
c     ipelet - 0 for no pellet; 1 to inject pellet
c     npel   - pellet counter
c     nupel  - I/O channel for pellet plot data (common /io/)
c ----------------------------------------------------------------------
c

c
      USE param
      USE fusion
      USE io 
      USE ions      
      USE neut
      USE nub
      USE nub2
      USE solcon
      USE soln
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE events
      USE pelcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray341.f,v 1.18 2007/10/18 00:03:05 stjohn Exp $"/
c
c      include 'param.i'
c      include 'events.i'
c      include 'fusion.i'
c      include 'geom.i'
c      include 'io.i'
c      include 'ions.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'neut.i'
c      include 'nub.i'
c      include 'nub2.i'
c      include 'numbrs.i'
c      include 'pelcom.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'sourc.i'
c
c     kj2 = 2*kj
c     kbg = max. number of fast beam ion energy groups
c     kag = max. number of fast alpha energy groups
c
      parameter (kj2 = kj*2, kbg = 25, kag = 2)
      dimension  dvol(kj),elden(kj),teev(kj),pden(kj),workp(kj)
      dimension  mapil(kj2),rchord(kj2)
      dimension  ebg(kbg),vcb(kj),denbig(kj,kbg),denbm(kbg)
      dimension  eag(kag),vca(kj),denaig(kj,kag)
c
c     check array dimensions
c
      if (ipelet .eq. 0)  return
      if (  2*nj .gt. kj2   )
     .  call STOP ('subroutine PROPEL: problem #1', 59)
      if (nbgpel .gt. kbg   )
     .  call STOP ('subroutine PROPEL: problem #2', 60)
      ibg = nbgpel/6
      if ( ibg*6 .ne. nbgpel)
     .  call STOP ('subroutine PROPEL: problem #3', 61)
c
c     set output unit numbers
c
      kprint = 3
      nprt   = nqik
      ndat   = nupel
c
c     increment pellet counter and set next pellet time
c
      npel = npel + 1
      write  (ndat, '(i8, e16.8)')  npel, time
      write  (ncrt,          90  )  npel, time
   90 format (' pellet number ', i3, ' at time = ', f7.3, ' seconds')
      timevent(8) = timpel(npel+1)
      if (npel .eq. 10)  timevent(8) = 1000.0
c
c          ONETWO has nj grid points enclosing njm1 plasma shells.
c          the pellet is assumed to pass horizontally through
c          the plasma center, thus it sees nseg = 2*njm1 plasma segments.
c          initialize plasma shell volumes, and map pellet path
c          segments onto the plasma shells.
c
      njm1 = nj-1
      nseg = 2*njm1
      pi   = 3.1415926
      do 100 j=1,njm1
 100  dvol(j)  = 2.0 * pi**2 * rmajor*(r(j+1)**2-r(j)**2)
      dvol(nj) = 0.0
      do 110 i=1,njm1
      mapil(nj-i)   = i
 110  mapil(njm1+i) = i
      do 120 i=1,nj
 120  rchord(i) = (r(nj)-r(nj+1-i)) / SQRT (kappa)
      do 130 i=nj+1,nseg+1
 130  rchord(i) = (r(nj)+r(i-nj+1)) / SQRT (kappa)
c
c          the pellet model uses the following strange grid:
c          te(1)  =  te(r=0)
c          te(i) = te(r = (r(i)+r(i+1))/2)
c          te(nj)  =   te(r=a)
c          convert ONETWO profiles to this grid.
c
      elden(1)  = ene(1)
      teev(1)   = te(1) * 1000.0
      do 140 i=2,njm1
      elden(i)  = 0.5*(ene(i)+ene(i+1))
 140  teev(i)   = 0.5*(te(i)+te(i+1)) * 1000.0
      elden(nj) = ene(nj)
      teev(nj)  = te(nj) * 1000.0
c
c          initialize beam fast ion energy groups.
c          there are nbgpel energy intervals of width debg from
c          the full injection energy down to zero.  nbgpel should
c          be a multiple of 6 to allow a clean match with the
c          half and third beam energy components.  note that
c          the interval from debg to 0 is always left empty.
c
      ebg(1)     = ebkev(1) * 1000.0
      debg       = ebg(1)/nbgpel
      nbg        = nbgpel-1
      do 150 jg=2,nbg
 150  ebg(jg)    = ebg(1) - (jg-1)*debg
      ebg(nbg+1) = debg
c
c     at each point in space, sum the beam particle density
c     from each energy component of each beamline.
c     model the fast ion slowing down by an analytical expression
c     derived by Callen and Rome, assuming steady state conditions.
c
      call zeroa (denbig,kj*kbg)
      call zeroa (denbm,kbg)
      do j=1,nj
      vcb(j) = 1.385e6 * SQRT (14.8*teev(j))
      do jbm=1,nbeams
      do 210 jbe=1,3
      nbge = nbgpel/jbe - 1
      jb1  = nbg - nbge + 1
      enbx = 0.5*(enb(j,jbe,jbm)+enb(j+1,jbe,jbm))
      if (j .eq. 1 .or. j .eq. nj)  enbx = enb(j,jbe,jbm)
      sum  = 0.0
      if (enbx .eq. 0.0)  go to 210
      do jbg=jb1,nbg
          rmass  = atw_beam*1.6726e-24
          rnfdot = sbsav(j,jbe,jbm)
          engy   = ebg(jbg)
****      e0     = ebg(1)/jbe
****      encrit = 0.5*atw_beam*1.6726e-24*1.6022e-12*vcb(j)**2
          v0     = SQRT (2.0*ebg(1)*1.6022e-12/(jbe*rmass))
          vgy    = SQRT (2.0*engy*1.6022e-12/rmass)
          tscat  = taus(j)
          svncx  = (enn(j,1)+enn(j,2))*cxr(engy*0.001/atw_beam)
          taucx  = 1.0e30
          if (svncx .gt. 0.0)  taucx = 1.0 / svncx
          denbm(jbg) = (rnfdot*tscat/rmass) * vgy *
     .         (1.0/(vgy**3+vcb(j)**3)) *
     .         (((v0**3+vcb(j)**3)/(vgy**3+vcb(j)**3))**(tscat/taucx)) *
     .         debg*1.6022e-12
          sum = sum + denbm(jbg)
      end do
c
      rksum = enbx/sum
      do jbg=jb1,nbg
        denbig(j,jbg) = denbig(j,jbg) + denbm(jbg)*rksum
      end do
  210 continue
      end do
      end do
c
c     neglect fast alpha particles for now
c
      nag    = 1
      atwalp = 4.0
      call zeroa (eag,kag)
      call zeroa (vca,kag)
      call zeroa (denaig,kag*kj)
c
c     fling the pellet
c
      call pellet (kprint, nprt, ndat, atw(ipel), pelrad, vpel,
     .             nj, kj, dvol, nseg, mapil, rchord, elden, teev,
     .             atw_beam, nbg, ebg, vcb, denbig,
     .             atwalp, nag, eag, vca, denaig, pden)
c
c          convert pden to ONETWO grid
c
      call copya (pden, workp, nj)
      pden(1)  = workp(1)
      pden(2)  = (workp(1) + 2.0*workp(2))/3.0
      do 250 i=3,njm1
 250  pden(i)  = 0.5 * (workp(i-1) + workp(i))
      pden(nj) = workp(nj)
c
c          calculate new plasma ene and te profiles.
c          note: 32.6 ev/pellet particle is dissipated in melting
c          and ionizing the pellet.  if fast beam ions are present,
c          it is assumed that they do all of the ablation, but the
c          ablation energy is not subtracted from the fast ion energy.
c
      eion = 0.0326
c
      do j=1,nj
        if (enb(1,1,1) .ne. 0.0)
     .  te(j)  = te(j)*ene(j)/(ene(j)+pden(j))
        if (enb(1,1,1) .eq. 0.0)
     .  te(j)  = (te(j)*ene(j)-eion*pden(j))/(ene(j)+pden(j))
        ene(j) = ene(j) + pden(j)
      end do
c
c     calculate new plasma ion density and temperature profiles
c
      do 310 j=1,nj
 310  workp(j) =  4.0 * pi**2 * rmajor*hcap(j)
      call trapv (r,en(1,ipel),workp,nj,totnum)
      do j=1,nj
        denion = 0.0
        do 320 i=1,nion
 320    denion     = denion + en(j,i)
        ti(j)      =  ti(j)*denion/(denion+pden(j))
        en(j,ipel) =  en(j,ipel) + pden(j)
      end do
c
      call trapv (r,en(1,ipel),workp,nj,totnew)
c
c     update the total number of ions for species ipel
c
      totadd = totnew - totnum
      if (ineut(ipel) .ne. 0)  snaddt(ipel) =  snaddt(ipel) + totadd
      return
c
      end

      real*8 function ssumpl (n, sx, incx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- ssumpl adds up the elements of the array sx.
c --- last revision: 3/81 w.a.houlberg and s.e.attenberger ornl.
c --- other comments:
c --- use OMNILIB version for CRAY optimization
c ----------------------------------------------------------------------
c
      dimension sx(*)
c
      sum = 0.0
      do i=1,n,incx
        sum = sum + sx(i)
      end do
      ssumpl = sum
      return
c
      end
