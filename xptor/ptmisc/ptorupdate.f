ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ptorupdate
c
c Derived from ptorinit_dv.f G.M.Staebler 7-8-00
c Updates geometry and Sources
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/ptor.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/vtrans.m'
      include '../inc/glf.m'
c
      real*8 volint  ! volume integrating function
      real*8 volintk ! break down grid volume int. func.
      real*8 vol(mxgrid,2)  !volume enclosed by zone center, zone boundary
      integer i,j,k
      integer icount
      integer kgrid
      real*8 powe_p,powe_mm,powi_p,powi_mm
      real*8 vbeam, cosbeam
      real*8 sour_p,sour_mm,mparsour_p,mparsour_mm
      real*8 lnlam,r_eps,fcm,kpol
      real*8 thetam,rhom,rminm
c
c Set default values for (some) namelist variables:
c
      chiip_mult=1.D0
      eigen_gf=0
c
c set a few useful constants
      pi=dabs(dacos(-1.D0))
      kJpereV=1.6022D-19   ! Boltzmann_s in Joules per eV.
c
c set weights
      do i=1,mxfields
        do k=1,ngrid
          vrho(i,k)=1.D0
        enddo
      enddo
c     
c setup up the geometry:
      do k=1,mxgrid
         r(k,2)=arho_exp*rho(k)
         r(k,1)=arho_exp*(rho(k)+rho(k-1))/2.D0
      enddo
c
      do k=1,mxgrid
         vol(k,1)=(vol_exp(k)+vol_exp(k-1))/2.D0
         vol(k,2)=vol_exp(k)
      enddo
c
c dr(k,1) = dr_k = width of zone k
c dr(k,2) = dr_{k+1/2} = distance between zone k and k+1
c
        do k=2,mxgrid-1
           dr(k,1)=r(k,2)-r(k-1,2)
           dr(k,2)=r(k+1,1)-r(k,1)
        enddo
        k=1
           dr(k,1)=r(1,2)
           dr(k,2)=r(k+1,1)-r(k,1)
        k=mxgrid
           dr(k,1)=r(k,2)-r(k-1,2)
c
      do k=2,mxgrid
         vprime(k,1)=(vol(k,2)-vol(k-1,2))/dr(k,1)
      enddo
      k=1
         vprime(k,1)=vol(k,2)/dr(k,1)
c
      do k=1,mxgrid-1
         vprime(k,2)=(vol(k+1,1)-vol(k,1))/dr(k,2)
      enddo
c      
c this is a fix to force ptor into agreement with modelt
c
      if(ifix.eq.1) then
       do k=1,mxgrid
        vprime(k,2)=vprime(k,1)
       enddo
      endif
      if(ifix.eq.2) then
       do k=1,mxgrid-1
        vprime(k,1)=vprime(k,2)
       enddo
      endif
      do k=1,mxgrid-1
       vphi_exp(k) = vphi_exp(k)/1.D3
       vpar_exp(k) = vpar_exp(k)/1.D3
       theta_exp(k)= arho_exp*rho(k)/(q_exp(k)*rmajor_exp)
         r_eps = (rmin_exp(k+1)+rmin_exp(k))
     >        /(rmaj_exp(k+1)+rmaj_exp(k))
         fcm = 1.D0+(0.46D0*r_eps-1.46D0)*DSQRT(r_eps)
         kpol=0.8839D0*fcm/(0.3477D0+0.4058D0*fcm)
cgms initialize vper_m to its neoclassical value (km/sec)
         vper_exp(k)=-alpha_dia/(ni_exp(k)*bt_exp*dr(k,2))*
     > ((1.D0-kpol)*(ni_exp(k+1)*ti_exp(k+1)-ni_exp(k)*ti_exp(k))
     > + kpol*(ti_exp(k+1)+ti_exp(k))/2.D0*(ni_exp(k+1)-ni_exp(k)))
     > - theta_exp(k)*vphi_exp(k)
      enddo
      theta_exp(0)=0.D0
      vper_exp(0)=vper_exp(1)
      vper_exp(mxgrid)=vper_exp(mxgrid-1)
      vphi_exp(mxgrid) = vphi_exp(mxgrid)/1.D3
      vpar_exp(mxgrid) = vpar_exp(mxgrid)/1.D3
      if (dvflag.eq.0) then
        do k=1,mxgrid
          vphi_exp(k) = vphi_exp(k)*1.D3
          vpar_exp(k) = vpar_exp(k)*1.D3
        enddo
      endif
c
c update the outer grid region
      do i=ngrid,mxgrid
        ni_m(i)=ni_exp(i)
        ne_m(i)=ne_exp(i)
        nz_m(i)=nz_exp(i)
        ti_m(i)=ti_exp(i)
        te_m(i)=te_exp(i)
        vphi_m(i)=vphi_exp(i)
        vper_m(i)=vper_exp(i)
        vpar_m(i)=vpar_exp(i)
      enddo
c update the non-evolved profiles
      if(itport_pt(1).eq.0)then
        do i=1,ngrid-1
          ni_m(i) = ni_exp(i)
          ne_m(i) = ne_exp(i)
          nz_m(i) = nz_exp(i)
          denpro(i)=0.1D0*ne_exp(i)
          dinpro(i)=0.1D0*ni_exp(i)
          dznpro(i)=0.1D0*nz_exp(i)
        enddo
        ni_m(0)=ni_m(1)
        ne_m(0)=ne_m(1)
        nz_m(0)=nz_m(1)
      endif
      if(itport_pt(2).eq.0)then
        do i=1,ngrid-1
          te_m(i) = te_exp(i)
        enddo
        te_m(0)=te_m(1)
      endif
      if(itport_pt(3).eq.0)then
        do i=1,ngrid-1
          ti_m(i) = ti_exp(i)
        enddo
        ti_m(0)=ti_m(1)
      endif
      if(itport_pt(4).eq.0)then
        do i=1,ngrid-1
          vphi_m(i) = vphi_exp(i)
          vpar_m(i) = vpar_exp(i)
          vphipro(i)=dabs(vphi_exp(i))
          vparpro(i)=vpar_exp(i)
        enddo
        vphi_m(0)=vphi_m(1)
      endif
      if(itport_pt(5).eq.0)then
        do i=1,ngrid-1
          vper_m(i) = vper_exp(i)
          vperpro(i)=vper_exp(i)
        enddo
        vper_m(0)=vper_m(1)
      endif
c
c Setups
c      
      Btor=bt_exp
      kappa=elonga_exp
      Rmaj=rmajor_exp
      amin=arho_exp/dsqrt(kappa)
      qa=q_exp(jmaxm)
      q0=q_exp(0) 
c compute ptor sources     
       do k=1,mxgrid-1
         zptheta_exp(k)=-2.D0*(theta_exp(k+1)-theta_exp(k))/
     >   (dr(k,2)*(theta_exp(k+1)+theta_exp(k)))
         powe_p=powe_exp(k)
         if(iexch.ge.1) powe_p=powe_p+pow_ei_exp(k)
         if(iohm.ge.1) powe_p=powe_p-powe_oh_exp(k)
         powe_mm=powe_exp(k-1)
         if(iexch.ge.1) powe_mm=powe_mm+pow_ei_exp(k-1)
         if(iohm.ge.1) powe_mm=powe_mm-powe_oh_exp(k-1)         
         powi_p=powi_exp(k)
         if(iexch.ge.1) powi_p=powi_p-pow_ei_exp(k) 
         powi_mm=powi_exp(k-1)
         if(iexch.ge.1) powi_mm=powi_mm-pow_ei_exp(k-1)
         Peaux(k)=(powe_p-powe_mm)/vprime(k,1)/dr(k,1)
         Piaux(k)=(powi_p-powi_mm)/vprime(k,1)/dr(k,1)
         Pohpro(k)=0.D0
         if(iohm.ge.1)Pohpro(k)=(powe_oh_exp(k)-powe_oh_exp(k-1))
     >      /vprime(k,1)/dr(k,1)
         Pe_alpha(k)=(powe_fus_exp(k)-powe_fus_exp(k-1))
     >      /vprime(k,1)/dr(k,1)
         Pi_alpha(k)=(powi_fus_exp(k)-powi_fus_exp(k-1))
     >      /vprime(k,1)/dr(k,1)
c
c pow_e,i_exp (total power) in MW= MJ/sec
c flow_exp in MW/kev=kamps 
         sour_p=flow_exp(k)
         sour_mm=flow_exp(k-1)
c Psour in MW/kev/m**3
         Psour(k)=(sour_p-sour_mm)/vprime(k,1)/dr(k,1) 
c Mpar and Mper momentum velocity sources
c cosbeam is a barm angle
c vbeam is a beam ion velocity in m/sec..same units as vpar, vper
c temporary: cosbeam=0.5   vbeam=beam velocity of 70kev beam
         cosbeam=0.5D0
         vbeam=1.D-2*9.79D5*
     >     DSQRT(2.D0*70.D0*1.D3/amassgas_exp)
         mparsour_p=flow_beam_exp(k)*cosbeam*vbeam
         mparsour_mm=flow_beam_exp(k-1)*cosbeam*vbeam
         Mpar(k)=pscale*cosbeam*vbeam*kevdsecpmw*sbeam_d(k+1)
         Mper(k)=0.D0
         nu_ei(k)=0.D0
      enddo
      zptheta_exp(0)=-2.D0*theta_exp(1)/(dr(1,2)*theta_exp(1))
      zptheta_exp(mxgrid)=zptheta_exp(mxgrid-1)
c
      do k=1,ngrid
       kgrid=k
       psume(k)=volintk(Peaux,kgrid)
       psumi(k)=volintk(Piaux,kgrid)
      enddo
c
crew6.30
      do j=0,jmaxm
       pow_ei_cor_m(j)=0.D0
      enddo
c
c jek - update boundary conditions
c
      denpro(ngrid)=0.1D0*ne_exp(ngrid)*(1.+didledge_n)
      dinpro(ngrid)=0.1D0*ni_exp(ngrid)*(1.+didledge_n)
      dznpro(ngrid)=0.1D0*nz_exp(ngrid)*(1.+didledge_n)
      Tipro(ngrid)=ti_exp(ngrid)*(1.+didledge_ti)
      Tepro(ngrid)=te_exp(ngrid)*(1.+didledge_te)
      vphipro(ngrid)=dabs(vphi_exp(ngrid))*(1.+didledge_vphi)
      vparpro(ngrid)=vpar_exp(ngrid)
      vperpro(ngrid)=vper_exp(ngrid)
c
      return
      end
