      subroutine palpha_dv
c
c same as palpha but uses dv method variables
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Calculate the alpha heating power using either revised
c NRL formula, Boucher PRETOR formula, or Bosch-Hale
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      implicit none
c
      include '../inc/ptor.m'
      include '../inc/tport.m'
c
      real*8 volint  ! volume integrating function
      integer k
      real*8 fusrate          ! fusion rate (/m**3/sec)
      real*8 Pfuspro(mxgrid)  ! fusion power density (MW/m**3)
      real*8 zTi,sigv,ecrit,frac,sigvv,sigvth
      real*8 denIO, nine, Tii, zalpha, massalpha, npne, 
     &       rtrit, cst1, theta, f, Wwcrit, dtrate_dv, dtrate_simple
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      if(ifusmodel.eq.1) then
c
c D. Boucher expression for <sigma v>_DT from PRETOR code
c
      zalpha = 2.D0              ! He charge
      massalpha = 4.D0           ! He mass
      rtrit = 0.5D0              ! nT/(nD+dT)
      cst1  = (1.D0-rtrit)*rtrit
      do k=1,ngrid
        denIO = ni_m(k)
        nine  = ni_m(k)/ne_m(k)
        Tii   = max(1.D-4,ti_m(k))
        npne  = nz_m(k)/ne_m(k)
        theta = Tii /(1.D0 - Tii * (1.51361D-2 + Tii *
     &          (4.60643D-3 - Tii * 1.0675D-4))/
     &          (1.D0 + Tii * ( 7.51886D-2 + Tii * (1.35D-2 +
     &          Tii * 1.366D-5))))
        f = (1182.176936D0/(4.D0 * theta))**(1.D0/3.D0)
        sigv = 1.17302D4 * 1.D-19 * theta *
     &         dsqrt(f/(1124656.D0 * Tii* Tii * Tii)) * 
     &         dexp(-3.D0 * f)
        Pfuspro(k)=fuscale*5.D0*5.6D0*cst1*sigv*1.D19*denIO**2.D0
        Wwcrit=58.3333D0/(te_m(k)*((nine/2.D0+
     &         zalpha*zalpha*npne/massalpha)**(2.D0/3.D0)))
        frac=1.03D0/(1.D0+Wwcrit*(0.35D0+0.02D0*Wwcrit))
        if(p_pulse_pt.eq.0.) then
         Pi_alpha(k)=      frac *Pfuspro(k)*3.5D0/17.5D0
         Pe_alpha(k)=(1.D0-frac)*Pfuspro(k)*3.5D0/17.5D0
        endif
c       write(*,50) k,rho(k),Tii,denIO,
c    &              sigv,Pfuspro(k)
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      elseif(ifusmodel.eq.2) then
c
c This simple expression for <sigma v>_DT is from the NRL plasma
c Formulary, which claims it is only good below 25 keV.  It is 
c about 30% low at 10 keV ...
c J. Kinsey 3/15/01 Constructed new fit that yields a
c better fit for temperatures between 7 and 100 keV. 
c     sigmav-DT=1.20e-13*T^-0.80*exp(-23.18*T^-.655)
c See comparison table below:
c
c    Ti(keV)    sigma-v(DT)     NRL fit     JEK fit 
c       7        4.07e-17      2.989e-17   3.880e-17
c      10        1.10e-16      7.580e-17   1.125e-16
c      20        4.20e-16      3.222e-16   4.202e-16
c      50        8.70e-16      1.209e-15   8.783e-16
c     100        8.50e-16      2.327e-15   9.686e-16
c
c From 7 kev down to 1 kev, a supplemental fit was constructed 
c yielding,
c     sigmav-DT=1.63e-13*T^0.20*exp(-17.21*T^-.352)
c
c    Ti(keV)    sigma-v(DT)     NRL fit     JEK fit 
c       1        5.50e-21      8.054e-21   5.470e-17
c       2        2.60e-19      3.104e-21   2.608e-19
c       5        1.30e-17      1.085e-27   1.289e-17
c       7        4.07e-17      2.989e-17   4.106e-17
c
c Setting idt=1 uses deuterium and tritium densities from
c ufiles if not a 50:50 mix.
c Note for NRL: 3.68e-12 in units cm3/s
c               so use 3.68e-18 for m3/s
c               and ni=ni_exp*1e19
c Total fusion rate: e19*e19*e-18=e20
c (the original factor e22 comes from dinpro=0.1*denpro*ni_exp/ne_exp)
c E_crit assuming pure 50:50 D:T plasma.  Assume alpha slows down
c only on electrons above E_crit, and only on ions below E_crit.
c
      do k=1,ngrid
         zTi=max(1.D-4,ti_m(k))
         if (zTi.le.7.) then
           sigv=1.63D-19*zTi**(0.20D0)*
     &          dexp(-17.21D0/(zTi**0.352D0))
         elseif (zTi.gt.7. .and. zTi.le.100.) then
           sigv=1.20D-19*zTi**(-0.80D0)*
     &          dexp(-23.18D0/(zTi**0.655D0))
         else
           sigv=3.68D-18/zTi**(2.D0/3.D0)*
     &          dexp(-19.94D0/zTi**(1.D0/3.D0))
         endif
c
c        sigv  = dtrate(zTi)*1.D-6   ! Bosch-Hale sigma-v
c
         if(idt.eq.1) then
           fusrate=(nm2_exp(k)*1.D19)*(nm3_exp(k)*1.D19)*sigv
         else
           fusrate=(ni_m(k)*1.D19/2.D0)**2.D0*sigv
         endif
         Pfuspro(k)=fuscale*fusrate*17.6D0*kJpereV
         ecrit=33.D0*te_m(k)
         frac=ecrit/3500.D0
         if(p_pulse_pt.eq.0.) then
           Pi_alpha(k)=      frac *Pfuspro(k)*3.5D0/17.5D0
           Pe_alpha(k)=(1.D0-frac)*Pfuspro(k)*3.5D0/17.5D0
         endif
c        write(*,50) k,rho(k),ti_m(k),ni_m(k)*1.D19,
c    >               sigv,Pfuspro(k),Pe_alpha(k)
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      elseif(ifusmodel.eq.3) then
c
c To match TGYRO
c
      do k=1,ngrid
         zTi=max(1.D-4,ti_m(k))
         sigv  = dtrate_simple(zTi)*1.D-6
         fusrate=(ni_m(k)*1.D19/2.D0)**2.D0*sigv
         Pfuspro(k)=fuscale*fusrate*17.5D0*kJpereV
         frac=0.2D0
         Pi_alpha(k)=frac*Pfuspro(k)*3.5D0/17.5D0
         Pe_alpha(k)=(1.D0-frac)*Pfuspro(k)*3.5D0/17.5D0
c         write(*,50) k,rho(k),ti_m(k),ni_m(k)*1.D19,
c     >               sigv,Pfuspro(k),Pe_alpha(k)
      enddo
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      else
c
c Use Bosch-Hale rate coefficient
c or TRANSP coefficients by using sigvv
c
      do k=1,ngrid
         zTi=max(1.D-4,ti_m(k))
         sigv  = dtrate_dv(zTi)*1.D-6   ! Bosch-Hale sigma-v
c         sigvv = sigvth(zTi)         ! TRANSP sigma-v
c          write(*,50) k, rho(k), ti_m(k), sigv
         if(idt.eq.1) then
           fusrate=(nm2_exp(k)*1.D19)*(nm3_exp(k)*1.D19)*sigv
         else
           fusrate=(ni_m(k)*1.D19/2.D0)**2.D0*sigv
         endif
         Pfuspro(k)=fuscale*fusrate*17.6D0*kJpereV
         ecrit=33.D0*te_m(k)
         frac=ecrit/3500.D0
         if(p_pulse_pt.eq.0.) then
           Pi_alpha(k)=      frac *Pfuspro(k)*3.5D0/17.5D0
           Pe_alpha(k)=(1.D0-frac)*Pfuspro(k)*3.5D0/17.5D0
         endif
c         write(*,50) k,rho(k),ti_m(k),ni_m(k)*1.D19,
c     >               sigv,Pfuspro(k),Pe_alpha(k)
      enddo
c
      endif
c
      Pfusion=volint(Pfuspro)
c
c Artificially prevent PFusion from exceeding PFusion_max.
c This gets around the thermonuclear instability.
c This limit could be achieved by feedback of auxiliary heating power,
c by feedback of the toroidal field ripple, or self-regulated by
c increasing transport losses as the temperature rises...
c
c Note that only the alpha heating power is limited, the actual
c total Pfusion will still be reported...
c
      if(Pfusion .gt. pfusion_max) then
         frac=pfusion_max/Pfusion
         do k=1,ngrid
            Pi_alpha(k)=frac*Pi_alpha(k)
            Pe_alpha(k)=frac*Pe_alpha(k)
         enddo
      endif
c
 50   format(i2,2x,0p1f5.3,0p6e13.5)
c
      return
      end
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function dtrate_dv (ti)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ------------------------------------------------- 11/20/95 --- HSJ ---
c     returns rate(cm**3/sec) of t(d,n)he4 reaction
c     New Bosch & Hale rate coefficient:
c     Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
c ----------------------------------------------------------------------
c
c     data for t(d,n)he4:
c
      data  C1,C2,C3,C4,C5,C6,C7, B_gsq, mrcsq
     .    / 1.17302D-9, 1.51361D-2, 7.51886D-2, 4.60643D-3,1.3500D-2,
     .     -1.06750D-4, 1.36600D-5,   1.182170D+3, 1124656.D0 /
      real*8 ti
c
c     neutrons produced by bulk plasma d-t fusion:
c
      theta  = ti*(C2 +ti*(C4+ti*C6))
      theta  = theta/(1.D0+ti*(C3+ti*(C5+ti*C7)))
      theta  = ti/(1.D0-theta)
      xsi    = (B_gsq/(4.D0*theta))**(1.D0/3.D0)
      dtrate_dv = C1*theta*DSQRT(xsi/(mrcsq*ti**3)) * DEXP(-3.D0*xsi)
      return
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      real*8 function dtrate_simple (ti)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
      real*8 r0,a1,a2,a3,a4,a5,a6
c
      r0 = 0.2935D0
      a1 = -21.378D0
      a2 = -25.204
      a3 = -7.101D-2
      a4 = 1.938D-4
      a5 = 4.925D-6
      a6 = -3.984D-8
c
      dtrate_simple = exp(a1/ti**r0+a2+ti*(a3+ti*(a4+ti*(a5+ti*a6))))
c
      return
c
      end
