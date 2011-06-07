      subroutine ip_chi1(itest_chi,RLT,RLN,q,shat,zth,nbeam,tau,eps,nu,
     .     rmajor,rho_i,v_ti,RLTcrit,RLTcrit2,chi_0,g,chi_i,chi_e)

c
c Definition of arguments:
c
c     RLT  R/L_Ti, where R = rmajor and L_Ti = -1/T_i dT_i/dr
c     RLN  R/L_ni, where L_ni = -1/n_i dn_i/dr
c     q    The usual.
c     shat == r/q dq/dr
c     zth  Thermal Z_eff.  The simulations that were carried out to 
c          generate the formulae in this subroutine assumed the plasma 
c          was composed of a thermal hydrogenic species, thermal carbon, 
c          a hydrogenic beam species, and electrons.  We found that high-Z
c          impurities primarily act to dilute the main ion concentration, 
c          and can accounted for to first order by modifying the 
c          definition of Z_eff.  Some of the more important effects of 
c          the fast ions in the plasma are also partially accounted 
c          for by this parameter:
c          zth == (n_i + 36 n_C)/(1 - n_beam)
c          Non-carbon light impurities could probably be substituted in 
c          the obvious way, but the threshold to the carbon-like branch
c          in the formulae below would then be only qualitatively accurate.  
c     nbeam == local fast ion (beam) density.
c     tau  == T_i/T_e
c     eps  == r/R
c     nu   Dimensionless collisionality parameter.
c          nu == 6.e-7 * n_e / (T_e**1.5 T_i**0.5) * rminor/80.
c          where n_e is in units of cm**-3, T_e and T_i are in eV, and 
c          rminor is in units of cm.  For an 80 cm, 100 eV, 1e13 plasma, 
c          nu=600.
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c update 5/9/95
c     gnu   Dimensionless collisionality parameter.
c          gnu == 2.5e-7 * n_e / (T_e**1.5 T_i**0.5) * rmajor
c          where n_e is in units of cm**-3, T_e and T_i are in eV, and 
c          rmajor is in units of m.  For an R = 2.4 m, 100 eV, 1e13 plasma, 
c          gnu=600.
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     rmajor == major radius of the plasma
c     rho_i == local thermal gyroradius of thermal hydrogenic species.  
c     v_ti == sqrt(T_i/m_i) where T_i and m_i are the local thermal 
c          hydrogenic temperature and average thermal hydrogenic mass.
c
c     Units: The only dimensional parameters in the inputs are the major
c     radius, rho_i, and v_t.  Their units should be consistent; the chi_s
c     that are returned will be in units of rho_i**2 v_ti / rmajor.
c
c OUTPUT:
c     RLTcrit: R/L_Tcrit for ITG mode
c     RLTcrit2: R/L_Tcrit for long-wavelength ITG mode
c     chi_0: normalized chi
c     g: L_Tc/L_T

c     This parameterization of chi is not complete for supershot parameters.
c     The problems are: 
c     (1) trapped electron effects, which often have a low growth rate (and 
c     thus probably a low chi) but which can alter the stability boundary 
c     significantly for low collisionality.  At high collisionality these modes
c     are largely irrelevant (so L-modes are probably unaffected).
c     (2) Gradients of Z_th_eff, which can strongly change chi and LT_crit.
c     This is probably important only for r/a<=0.3 and r/a>=0.7 according to 
c     estimates based on discussion with Synakowski.  However, they can be 
c     especially significant in pellet shots.
c     (3) Sheared flows, which are large enough to stabilize the ITG modes
c     at around 0.25<r/a<0.35 for several shots studied so far.  Sheared flows
c     affect the critical gradient much more than chi0.(?)
c     
      implicit none

      integer itest_chi

      real*8 RLT,RLN,shat,zth,tau,rminor,rmajor,radius,chi_i,q,
     .     nbeam,taub,rho_i,v_ti,eps,chi_e,nu,chi_0,g
      real*8 RLTcrit,chi0,alpha,RLTcrit2,RLN_th
      real*8 c16,c17,c18,f_0
      real*8 g_fac1,g_fac2
      data alpha/1.D0/

c      common/ip_common/ alpha
      if(zth.lt.3.3) then
         c16=1.54D0-.47D0*(3.3D0-zth)**0.6D0
      else
         c16=.55D0+.99D0*dexp(20.D0*(3.3D0-zth))
      endif

      if(zth.le.3) c17=1.D0+0.1D0/tau
      if(zth.gt.3.0 .and. zth.le.3.5) c17=2.D0*(3.5D0-zth)+0.1D0/tau
      if(zth.gt.3.5) c17=(0.1D0+0.2D0*(zth-3.5D0))/tau

      c18=0.61D0-0.27D0/(1.D0+dexp(4.D0*(3.5D0-zth)))
c
      taub=tau/(1.D0-nbeam)
c
C     critical ion temperature gradient:
C
      RLTcrit=2.778D0*dsqrt(0.5D0+1.D0/q)*c16*taub**c18
     .     *(abs(0.1976D0-0.4550D0*shat+0.1616D0*RLN)**0.769D0
     .     +0.7813D0+0.2762D0*shat+0.3967D0*shat**2)
     .     *abs(1.D0-0.85D0*eps/shat**0.25D0)
     .     *(1.D0-0.9D0*taub**0.4D0/nu**0.27D0*eps)
c rew note last line ....funny param at low nu???
c insert RLTcrit > 0 condition
      if(RLTcrit.le.0.) RLTcrit=1.D-6
C     
C     second critical ion temperature gradient:
C     
      RLN_th=max(0.D0,(RLn-8.7D0)/8.7D0)
      RLTcrit2=9.72D0*dsqrt(zth)/tau**0.4D0*(1.D0+RLN_th)
c          
      f_0=10.D0*q/(2.D0+shat)*c17/taub*(1.D0+6.D0*eps/nu**0.14D0)
      chi0=f_0 * rho_i**2*v_ti/rmajor     
c     
      if(RLT-RLTcrit.gt.0.) then
         g_fac1=min(dsqrt(RLT-RLTcrit),(RLT-RLTcrit))
      else
         g_fac1=0.D0
      endif
      
      if(alpha.eq.-1) then
         write(*,*) 'alpha=?'
         read(*,*) alpha
      endif
      if(RLT-RLTcrit2.gt.0.) then
         g_fac2=alpha*(RLT-RLTcrit2)
      else
         g_fac2=0.D0
      endif
      
      g=RLT/RLTcrit
      chi_0=chi0*dsqrt(RLTcrit)
      chi_i=chi0*(g_fac1+g_fac2)
      chi_e=chi0*g_fac1*0.22D0*tau**0.6D0*
     .      dsqrt(q/(1.5D0*shat))*nu**0.14D0
  
      if(alpha.eq.1.e-10) then
       write(*,*) q,shat,eps,nu,c17,taub
       write(*,*) f_0, chi0, RLT, RLTcrit, g_fac1, chi_0, chi_i
      endif
      
      return
      end
      



