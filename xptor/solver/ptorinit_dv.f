ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ptorinit_dv
c
c Read input data and initialize geometry and other variables
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
      real*8 vol(mxgrd,2)  !volume enclosed by zone center, zone boundary
      integer i,j,k
      integer icount
      integer kgrid,nsum
      real*8 powe_p,powe_mm,powi_p,powi_mm
      real*8 vbeam, cosbeam, cnc
      real*8 sour_p,sour_mm,mparsour_p,mparsour_mm
      real*8 lnu_eilam,r_eps,fcm,kpol1, kpol2, b_unit
      real*8 rhom,rminm,qm,temm,gradrm
      real*8 thetamin,Peaux_tot,Piaux_tot,Pohm_tot
      real*8 ugradte,ugradti,ugradne
      real*8 ugradfi,ugradfz,ugradni,ugradnz
      real*8 thetam,dvoldr
      real*8 gradr_exp(mxgrd)
      real*8 ctorm,grad_c_tor,gradvphim
      real*8 nuei,zbrac,tew,tiw,new,niw,nzw,lnlam
      real*8 vthe,vthi,vthz
c
c... GLF models
c    0 = original GLF23
c    1 = retuned GLF23 (v1.61)
c
      chiip_mult=1.D0
      eigen_gf=0
c
        if (iglf.eq.-1) then     ! renormalization only
          cnorm_gf=27.D0         ! ITG/TEM normalization
          xparam_gf(10)=17.8D0   ! ETG normalization (cnorm*xparam(10))
          xparam_gf(13)=0.2      ! rms_theta q-dependence
          xparam_gf(15)=-0.1     ! trapped ptcle fraction reduction
          xparam_gf(16)=0.0      ! rms_theta shat dependence
          xparam_gf(17)=0.10     ! rms_theta shat dependence
          xparam_gf(19)=0.0      ! rms_theta alpha-dependence
          iflagin_gf(5)=3        ! rms_theta
          park_gf=0.7D0
          adamp_gf=0.5D0
          alpha_p_gf=0.50D0
          bt_flag=0
        endif
        if (iglf.eq.1) then      ! v1.61 (retuned model)
          cnorm_gf=50.D0         ! ITG/TEM normalization
          xparam_gf(10)=12.D0    ! ETG normalization (cnorm*xparam(10))
          xparam_gf(13)=0.2      ! rms_theta q-dependence
          xparam_gf(15)=-0.1     ! trapped ptcle fraction reduction
          xparam_gf(16)=0.15     ! rms_theta shat dependence
          xparam_gf(17)=0.25     ! rms_theta shat dependence
          xparam_gf(19)=1.0      ! rms_theta alpha-dependence
          iflagin_gf(5)=6        ! rms_theta w/ tau**0.25
          park_gf=0.8
          adamp_gf=0.7
          alpha_p_gf=0.35
          bt_flag=1
        endif
c
c
c set a few useful constants
      pi=dabs(dacos(-1.D0))
      kJpereV=1.6022D-19   ! Boltzmann_s in Joules per eV.
      kevdsecpmw=1.6022D-19*1.D3*1.D-6
c
c set the weights of the vpol equation
      if(itport_pt(5).le.1)then
       wp(1) = 1.D0
       wp(2) = 0.D0
       wp(3) = 0.D0
      elseif(itport_pt(5).eq.2)then
       wp(1) = 0.D0
       wp(2) = 1.D0
       wp(3) = 0.D0
      elseif(itport_pt(5).eq.3)then
       wp(1) = 0.D0
       wp(2) = 0.5D0
       wp(3) = 0.5D0
      endif
c
c     
c setup up the geometry:
c         drhodr(0)=1.D0
c         geofac(0)=1.D0
c         drhodrrrho(0)=1.D0
c         gradrhosq_exp(0)=1.D0
      do k=1,mxgrid
         r(k,2)=arho_exp*rho(k)
         r(k,1)=arho_exp*(rho(k)+rho(k-1))/2.D0
c         drhodr(k)=1.D0
c         geofac(k)=1.D0
c         drhodrrrho(k)=1.D0
c         gradrhosq_exp(k)=1.D0
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
c
c Setups
c      
      Btor=bt_exp
      kappa=elonga_exp
      Rmaj=rmajor_exp
      amin=arho_exp/dsqrt(kappa)
      qa=q_exp(jmaxm)
      q0=q_exp(0) 
      cnc =  -1.0*alpha_dia/bt_exp 
      cv=1000.0    
c
c Initialize profiles:
c
       if(imodel.eq.82)call tglf_startup
c
c      thetamin=0.05
c... ion density profiles via dilution_model:
c    -1 : use ni_exp, nz_exp as set in datmap.f
c     0 : single impurity and fast ion dilution (default)
c     1 : fast ion dilution only
c     2 : ni=ne, nz=0
c
      do k=1,mxgrid
         nz_exp(k)=ne_exp(k)*(zeff_exp(k)-1.0D0)/
     >    (zimp_gf**2-zimp_gf)
         if(dilution_model.eq.0)then
c full imputity and fast ion dilution
           ni_exp(k)=ne_exp(k)-nfst_exp(k)-zimp_gf*nz_exp(k) 
         elseif(dilution_model.eq.1)then
c just fast ion dilution
           ni_exp(k) = ne_exp(k)-nfst_exp(k)
           nz_exp(k) = 0.0
         elseif(dilution_model.eq.2)then
c no dilution
           ni_exp(k) = ne_exp(k)
           nz_exp(k) = 0.0
         endif
         fi_m(k) = ni_exp(k)/ne_exp(k)
         fz_m(k) = nz_exp(k)/ne_exp(k)
         mask_r(k) = 1
         theta_exp(k)= arho_exp*rho(k)/(q_exp(k)*rmajor_exp)
c         if(theta_exp(k).lt.thetamin)theta_exp(k)=thetamin
c         bp_exp(k) = (psir_exp(k)-psir_exp(k-1))/dr(k,1)
c         bp_exp(k) = bt_exp*drhodrrrho(k)
         bp_exp(k) = bt_exp*theta_exp(k)
        if(imodel.ne.82)then
          xb2_exp(k) = (bt_exp**2)/(f_exp(k)*h_exp(k))
     >     + (bp_exp(k)**2)*g_exp(k)
c         xb2_exp(k) = 1.0/(h_exp(k)*f_exp(k))+
c     >    g_exp(k)*theta_exp(k)**2
c temporary: should be <R**2>
          xr2_exp(k) = rmaj_exp(k)**2
     >     +0.5D0*rmin_exp(k)**2
          c_par(k) = xb2_exp(k)/bt_exp**2
          c_tor(k) = xr2_exp(k)/rmajor_exp**2
          c_per(k) = 1.0/ABS(f_exp(k))
          a_pol(k) = 1.0
          a_tor(k) = 1.0
        endif
c         
c  convert velocities to km/s
c
c  remember angrot_exp is the impurity ion rotation
c
         vphi_exp(k)=c_tor(k)*rmajor_exp*angrot_exp(k)/cv
c  
         nuei_m(k) = 0.0
         flow_m(k) = 0.0
         powe_m(k) = 0.0
         powi_m(k) = 0.0
         stress_tor_m(k) = 0.0
         stress_par_m(k) = 0.0
         pow_ei_cor_m(k) = 0.0
         stress_par_cor_m(k) = 0.0
         flow_glf(k) = 0.0
         powe_glf(k) = 0.0
         powi_glf(k) = 0.0
         stress_tor_glf(k) = 0.0
         stress_par_glf(k) = 0.0
         flow_neo(k) = 0.0
         powe_neo(k) = 0.0
         powi_neo(k) = 0.0
         stress_tor_neo(k) = 0.0
         stress_par_neo(k) = 0.0
         flow_adhoc(k) = 0.0
         powe_adhoc(k) = 0.0
         powi_adhoc(k) = 0.0
         stress_tor_adhoc(k) = 0.0
         stress_par_adhoc(k) = 0.0
         vpol_exp(k) = 0.0
         Pradb(k) = 0.0
         Prads(k) = 0.0
         Pi_alpha(k)=0.0
         Pe_alpha(k)=0.0
      enddo
      bp_exp(0)=0.0
      theta_exp(0)=0.0
      te_exp(0)=te_exp(1)
      ti_exp(0)=ti_exp(1)
      ne_exp(0)=ne_exp(1)
      ni_exp(0)=ni_exp(1)
      nz_exp(0)=nz_exp(1)
      vphi_exp(0)=vphi_exp(1)
      fi_m(0)=fi_m(1)
      fz_m(0)=fz_m(1)
      nuei_m(0) = nuei_m(1)
      flow_m(0) = 0.0
      powe_m(0) = 0.0
      powi_m(0) = 0.0
      stress_tor_m(0) = 0.0
      stress_par_m(0) = 0.0
      flow_glf(0) = 0.0
      powe_glf(0) = 0.0
      powi_glf(0) = 0.0
      stress_tor_glf(0) = 0.0
      stress_par_glf(0) = 0.0
      flow_neo(0) = 0.0
      powe_neo(0) = 0.0
      powi_neo(0) = 0.0
      stress_par_neo(0) = 0.0
      stress_tor_neo(0) = 0.0
      flow_adhoc(0) = 0.0
      powe_adhoc(0) = 0.0
      powi_adhoc(0) = 0.0
      stress_par_adhoc(0) = 0.0
      stress_tor_adhoc(0) = 0.0
      pow_ei_cor_m(0) = 0.0
      stress_par_cor_m(0) = 0.0
      vpol_exp(0) = 0.0
      c_par(0) = c_par(1)
      c_per(0) = c_per(1)
      c_tor(0) = c_tor(1)
      a_pol(0) = a_pol(1)
      a_tor(0) = a_tor(1)
c
c compute vexb_exp, vpar_exp and zptheta_exp
c
      a_unit_exp = rmin_exp(mxgrid)
c      if(igeo_tg.eq.0)a_unit_exp=arho_exp
c
      pow_ei_exp(0) = 0.0
      do k=1,mxgrid-1
        zptheta_exp(k)=-2.D0*(theta_exp(k+1)-theta_exp(k))/
     >   (dr(k,2)*(theta_exp(k+1)+theta_exp(k)))
        nem = (ne_exp(k+1)+ne_exp(k))/2.D0
        tim = (ti_exp(k+1)+ti_exp(k))/2.D0
        tem = (te_exp(k+1)+te_exp(k))/2.D0
        fim = (fi_m(k+1)+fi_m(k))/2.D0
        fzm = (fz_m(k+1)+fz_m(k))/2.D0
        nim = fim*nem
        nzm = fzm*nem
        vexbm= 0.0
        vpolm= 0.0
        gradnem = (ne_exp(k+1)-ne_exp(k))/dr(k,2)
        gradtim = (ti_exp(k+1)-ti_exp(k))/dr(k,2)
        gradtem = (te_exp(k+1)-te_exp(k))/dr(k,2)
        gradvexbm = 0.0
        gradvpolm = 0.0
        gradfim = (fi_m(k+1)-fi_m(k))/dr(k,2)
        gradfzm = (fz_m(k+1)-fz_m(k))/dr(k,2)
        gradnim = fim*gradnem + nem*gradfim
        gradnzm = fzm*gradnem + nem*gradfzm
        zpne_exp(k) = -a_unit_exp*(gradnem/nem)
        zpni_exp(k) = -a_unit_exp*(gradnim/nim)
        zpnz_exp(k) = -a_unit_exp*(gradnzm/nzm)
        zpte_exp(k) = -a_unit_exp*(gradtem/tem)
        zpti_exp(k) = -a_unit_exp*(gradtim/tim)
        jm=k
        call neoclassical
        do i=1,nspecies
          vneo_exp(i,k+1) = vneo(i)
          vdia_exp(i,k+1) = vdia(i)
        enddo
c compute energy exchange
        tew = te_exp(k)
        new = ne_exp(k)
        niw = ni_exp(k)
        nzw = nz_exp(k)
        tiw = ti_exp(k)
        zbrac=(niw+amassgas_exp*nzw*zimp_exp**2/amassimp_exp)/new
        lnlam=24.D0-DLOG(DSQRT(1.D13*new)/(1000.D0*tew))
        lnlam=DMAX1(lnlam,1.D0)
        nuei= 1.5D0*new*new*zbrac
     &     *3.2D-9*1.D13*lnlam
     &     /amassgas_exp/DSQRT(1000.D0*tew)**3
        pow_ei_exp(k) = pow_ei_exp(k-1) 
     >  -vprime(k,1)*dr(k,1)*1.6022D-3*nuei*(tiw-tew)
        if(iexch.eq.0)pow_ei_exp(k)=0.0
c
c compute power balance chi's
        dvoldr = vprime(k,2)*drhodr(k)*drhodr(k)
        diff_exp(k)= -gradnem*flow_exp(k)
     >  /(dvoldr*1.6022D-3*MAX(1.0D-12,gradnem**2))
        chie_exp(k)= -gradtem*(powe_exp(k)-pow_ei_exp(k))
     >   /(dvoldr*nem*1.6022D-3*MAX(1.0D-12,gradtem**2))
        chii_exp(k)= -gradtim*(powi_exp(k)+pow_ei_exp(k))
     >   /(dvoldr*nim*1.6022D-3*MAX(1.0D-12,gradtim**2))
      enddo
      diff_exp(0) = diff_exp(1)
      chie_exp(0) = chie_exp(1)
      chii_exp(0) = chii_exp(1)
      diff_exp(mxgrid) = diff_exp(mxgrid-1)
      chie_exp(mxgrid) = chie_exp(mxgrid-1)
      chii_exp(mxgrid) = chii_exp(mxgrid-1)
      do i=1,nspecies
        vdia_exp(i,1) = vdia_exp(i,2)
        vneo_exp(i,1) = vneo_exp(i,2)
        vdia_exp(i,0) = vdia_exp(i,1)
        vneo_exp(i,0) = vneo_exp(i,1)
      enddo
      k=mxgrid
        tew = te_exp(k)
        new = ne_exp(k)
        niw = ni_exp(k)
        nzw = nz_exp(k)
        tiw = ti_exp(k)
        zbrac=(niw+amassgas_exp*nzw*zimp_exp**2/amassimp_exp)/new
        lnlam=24.D0-DLOG(DSQRT(1.D13*new)/(1000.D0*tew))
        lnlam=DMAX1(lnlam,1.D0)
        nuei= 1.5D0*new*new*zbrac
     &     *3.2D-9*1.D13*lnlam
     &     /amassgas_exp/DSQRT(1000.D0*tew)**3
        pow_ei_exp(k) = pow_ei_exp(k-1) 
     >  -vprime(k,1)*dr(k,1)*1.6022D-3*nuei*(tiw-tew)
        if(iexch.eq.0)pow_ei_exp(k)=0.0
c
      do k=1,mxgrid-1
c initialize vexb_exp to its neoclassical value (km/sec)
c assuming ion species 3 is the CER measured one
         vexb_exp(k) = (vphi_exp(k)-vneo_exp(3,k)*c_per(k))/c_tor(k)
     >   - vdia_exp(3,k)
c reset vphi_exp to the main ion toroidal flow
c         vphi_exp(k) = c_per(k)*(vpol_exp(k)+vneo_exp(2,k))
c     >    + c_tor(k)*(vexb_exp(k)+vdia_exp(2,k))
c
c compute the parallel flows for each species
         vpar_exp(k)=a_pol(k)*vneo_exp(2,k)
     >    +a_tor(k)*(vdia_exp(2,k)+vexb_exp(k))
         tem = (te_exp(k+1)+te_exp(k))/2.0
         tim = (ti_exp(k+1)+ti_exp(k))/2.0
c thermal velocities in km/sec
         vthi = 9.79D3*DSQRT(2.D3*tim/amassgas_exp)/cv
         vthz = DSQRT(amassgas_exp/amassimp_exp)*vthi
         vthe = DSQRT(tem*amassgas_exp*1.8362D3/tim)*vthi
         mach_exp(1,k)=(a_pol(k)*vneo_exp(1,k)
     >    +a_tor(k)*(vdia_exp(1,k)+vexb_exp(k)))/vthe
         mach_exp(2,k)=(a_pol(k)*vneo_exp(2,k)
     >    +a_tor(k)*(vdia_exp(2,k)+vexb_exp(k)))/vthi
         mach_exp(3,k)=(a_pol(k)*vneo_exp(3,k)
     >    +a_tor(k)*(vdia_exp(3,k)+vexb_exp(k)))/vthz
       enddo
      vpar_exp(0) = vpar_exp(1)
c      vpar_exp(1,0) = vpar_exp(1,1)
c      vpar_exp(2,0) = vpar_exp(2,1)
c      vpar_exp(3,0) = vpar_exp(3,1)
      vexb_exp(0) = vexb_exp(1)
      vphi_exp(0) = vphi_exp(1)
      zptheta_exp(0)=zptheta_exp(1)
      vexb_exp(mxgrid) = vexb_exp(mxgrid-1)
      zptheta_exp(mxgrid)=zptheta_exp(mxgrid-1)
      vphi_exp(mxgrid) = vphi_exp(mxgrid-1)
      vpar_exp(mxgrid) = vpar_exp(mxgrid-1)
c
c  overwrites for restart
c
        do k=ngrid,mxgrid
         ne_m(k)=ne_exp(k)
         ni_m(k)=ni_exp(k)
         nz_m(k)=nz_exp(k)
         ti_m(k)=ti_exp(k)
         te_m(k)=te_exp(k)
         vexb_m(k)=vexb_exp(k)
         vpol_m(k)=0.0
         vphi_m(k)=vphi_exp(k)
         vpar_m(k)=vpar_exp(k)
         do i=1,nspecies
c           vphi_m(i,k)=vphi_exp(i,k)
c           vpar_m(i,k)=vpar_exp(i,k)
           vneo_m(i,k)=vneo_exp(i,k)
           vdia_m(i,k)=vdia_exp(i,k)
           mach_m(i,k)=mach_exp(i,k)
         enddo
         zpne_m(k) = zpne_exp(k)
         zpni_m(k) = zpni_exp(k)
         zpnz_m(k) = zpnz_exp(k)
         zpte_m(k) = zpte_exp(k)
         zpti_m(k) = zpti_exp(k)
        enddo
      if(restart_pt.le.0.0)then
        do k=1,ngrid
         ne_m(k)=ne_exp(k)*(1.0 + didledge_n)
          ni_m(k) = fi_m(k)*ne_m(k)
          nz_m(k) = fz_m(k)*ne_m(k)
         ti_m(k)=ti_exp(k)*(1.0 + didledge_ti)
         te_m(k)=te_exp(k)*(1.0 + didledge_te)
         vexb_m(k)=vexb_exp(k)
         vpol_m(k)=0.0
         vphi_m(k)=vphi_exp(k)
         vpar_m(k)=vpar_exp(k)
         vphi_m(k)=vphi_exp(k)
         vpar_m(k)=vpar_exp(k)
         do i=1,nspecies
c           vphi_m(i,k)=vphi_exp(i,k)
c           vpar_m(i,k)=vpar_exp(i,k)
           vneo_m(i,k)=vneo_exp(i,k)
           vdia_m(i,k)=vdia_exp(i,k)
         enddo
        enddo
      endif
c
      if(restart_pt.ge.0.0) then
        open (9,file='nepulse.dat',status='unknown')
        open (10,file='tepulse.dat',status='unknown')
        open (12,file='tipulse.dat',status='unknown')
        open (15,file='vexbpulse.dat',status='unknown')
        open (31,file='vpolpulse.dat',status='unknown')
c cgms       write(6,*)"restart_pt = ", restart_pt
        read(9,150)  ntime_t, ngrid
        read(10,150) ntime_t, ngrid
        read(12,150) ntime_t, ngrid
        read(15,150) ntime_t, ngrid
        read(31,150) ntime_t, ngrid
        do j=0,ntime_t
         read(9,310) time_t(j), ne_t(0:ngrid,j)
         read(10,310) time_t(j), te_t(0:ngrid,j)
         read(12,310) time_t(j), ti_t(0:ngrid,j)
         read(15,310) time_t(j), vexb_t(0:ngrid,j)
         read(31,310) time_t(j), vpol_t(0:ngrid,j)
         if(time_t(j).ge.restart_pt) exit
        enddo
        if(j.gt.ntime_t)j=ntime_t
        restart_pt=time_t(j)
        if(i_proc.eq.0)write(*,*)"restart at ",restart_pt,j
        do k=0,ngrid
          ne_m(k)=ne_t(k,j)
          ni_m(k)=fi_m(k)*ne_m(k)
          nz_m(k)=fz_m(k)*ne_m(k)
          te_m(k)=te_t(k,j)
          ti_m(k)=ti_t(k,j)
          vexb_m(k)=vexb_t(k,j)
          vpol_m(k)=vpol_t(k,j)
        enddo
        close(9)
        close(10)
        close(12)
        close(15)
        close(31) 
c
c compute vdia_m and vneo_m
c
       do k=1,ngrid-1
        nem = (ne_m(k+1)+ne_m(k))/2.D0
        tim = (ti_m(k+1)+ti_m(k))/2.D0
        tem = (te_m(k+1)+te_m(k))/2.D0
        fim = (fi_m(k+1)+fi_m(k))/2.D0
        fzm = (fz_m(k+1)+fz_m(k))/2.D0
        nim = fim*nem
        nzm = fzm*nem
        vexbm = 0.0
        vpolm = 0.0
        gradnem = (ne_m(k+1)-ne_m(k))/dr(k,2)
        gradtim = (ti_m(k+1)-ti_m(k))/dr(k,2)
        gradtem = (te_m(k+1)-te_m(k))/dr(k,2)
        gradvexbm = 0.0
        gradvpolm = 0.0
        gradfim = (fi_m(k+1)-fi_m(k))/dr(k,2)
        gradfzm = (fz_m(k+1)-fz_m(k))/dr(k,2)
        gradnim = fim*gradnem + nem*gradfim
        gradnzm = fzm*gradnem + nem*gradfzm
        zpne_m(k) = -a_unit_exp*gradnem/nem
        zpni_m(k) = -a_unit_exp*gradnim/nim
        zpnz_m(k) = -a_unit_exp*gradnzm/nzm
        zpte_m(k) = -a_unit_exp*gradtem/tem
        zpti_m(k) = -a_unit_exp*gradtim/tim
        jm = k
        call neoclassical
c
        do i=1,nspecies
          vneo_m(i,k+1) = vneo(i)
          vdia_m(i,k+1) = vdia(i)
        enddo
      enddo
c
      do k=1,ngrid
c compute the parallel flows for each species
         vpar_m(k)=a_pol(k)*(vneo_m(2,k)+vpol_m(k))
     >    +a_tor(k)*(vdia_m(2,k)+vexb_m(k))
c         vpar_m(1,k)=c_par(k)*(vneo_m(1,k)+vpol_m(k))
c     >    +c_per(k)*(vdia_m(1,k)+vexb_m(k))
c         vpar_m(2,k)=c_par(k)*(vneo_m(2,k)+vpol_m(k))
c     >    +c_per(k)*(vdia_m(2,k)+vexb_exp(k))
c         vpar_m(3,k)=c_par(k)*(vneo_m(3,k)+vpol_m(k))
c     >    +c_per(k)*(vdia_m(3,k)+vexb_m(k))
c  toridal imurity ion flow
         vphi_m(k) = c_per(k)*(vpol_m(k)+vneo_m(3,k))
     >    + c_tor(k)*(vexb_m(k)+vdia_m(3,k))
c
       enddo
      endif 
c set central derivatives to zero
      ne_m(0)=ne_m(1)
      te_m(0)=te_m(1)
      ti_m(0)=ti_m(1)
      vpol_m(0)=vpol_m(1)
      vexb_m(0)=vexb_m(1)
      vphi_m(0)=vphi_m(1)
      vpar_m(0)=vpar_m(1)
      do i=1,nspecies
        vneo_m(i,1) = vneo_m(i,2)
        vdia_m(i,1) = vdia_m(i,2)
c        vphi_m(i,0)=vphi_m(i,1)
c        vpar_m(i,0)=vpar_m(i,1)
        vneo_m(i,0)=vneo_m(i,1)
        vdia_m(i,0)=vdia_m(i,1)
      enddo
      ni_m(0)=ni_m(1)
      nz_m(0)=nz_m(1)
c 
c  boundary conditions
c
      vexb_m(ngrid)=vexb_m(ngrid) + vexb_bc
      vpol_m(ngrid)=0.0
c
c compute the line-integrated experimental electron density for feedback
c excluding sawtooth zone q < 1
c      wall_mult=1.0
c      te_bc_mult=1.0
c      ti_bc_mult=1.0
      sum_ne_exp=0.0
      sum_te_exp=0.0
      sum_ti_exp=0.0
      nsum=0
      do k=0,ngrid-1
        if(q_exp(k).gt.1.0)then
          sum_ne_exp = sum_ne_exp + ne_exp(k)
          sum_te_exp = sum_te_exp + te_exp(k)
          sum_ti_exp = sum_ti_exp + ti_exp(k)
          nsum = nsum + 1
        endif
      enddo
      sum_ne_exp = sum_ne_exp/REAL(nsum)
      sum_te_exp = sum_te_exp/REAL(nsum)
      sum_ti_exp = sum_ti_exp/REAL(nsum)
c
c  compute fusion, radiation ohmic and auxiliary powers
      if(ialpha.eq.1)call palpha_dv
      if(irad.eq.1)call prad_dv
      if(iohm.ge.1)call pohmic_dv
c
       flow_exp(0) = 0.0
       stress_tor_exp(0) = 0.0
       stress_par_exp(0) = 0.0
       do k=1,mxgrid-1
         powe_p=powe_exp(k)
c         if(iexch.ge.1) powe_p=powe_p+pow_ei_exp(k)
c         if(iohm.ge.1) powe_p=powe_p-powe_oh_exp(k)
         powe_mm=powe_exp(k-1)
c         if(iexch.ge.1) powe_mm=powe_mm+pow_ei_exp(k-1)
c         if(iohm.ge.1) powe_mm=powe_mm-powe_oh_exp(k-1)         
         powi_p=powi_exp(k)
c         if(iexch.ge.1) powi_p=powi_p-pow_ei_exp(k) 
         powi_mm=powi_exp(k-1)
c         if(iexch.ge.1) powi_mm=powi_mm-pow_ei_exp(k-1)
         Peaux(k)=(powe_p-powe_mm)/vprime(k,1)/dr(k,1)
         Piaux(k)=(powi_p-powi_mm)/vprime(k,1)/dr(k,1)
         Pohpro(k)=0.D0
c this is inconsistent?
c         if(iohm.ge.1)Pohpro(k)=(powe_oh_exp(k)-powe_oh_exp(k-1))
c     >      /vprime(k,1)/dr(k,1)
         if(ialpha.eq.0)then
           Pe_alpha(k)=(powe_fus_exp(k)-powe_fus_exp(k-1))
     >      /vprime(k,1)/dr(k,1)
           Pi_alpha(k)=(powi_fus_exp(k)-powi_fus_exp(k-1))
     >      /vprime(k,1)/dr(k,1)
           Pe_alpha(k)=0.0
           Pi_alpha(k)=0.0
         endif
c
c pow_e,i_exp (total power) in MW= MJ/sec
c flow_exp in KA 
         sour_p=(flow_exp(k)-flow_wall_exp(k))
         sour_mm=(flow_exp(k-1)-flow_wall_exp(k-1))
c Psour in KA/m^3
         Psour(k)=(sour_p-sour_mm)/(vprime(k,1)*dr(k,1))
c wall source
         sour_p=flow_wall_exp(k)
         sour_mm=flow_wall_exp(k-1)
         Psour_wall(k)=(sour_p-sour_mm)/(vprime(k,1)*dr(k,1))
c Mphi and Mpar momentum velocity sources
c cosbeam is a barm angle
c vbeam is a beam ion velocity in m/sec..same units as vpar
         cosbeam=0.5D0
c effective beam velocity from 0.5*m*vbeam^2 = power/flux
          vbeam = SQRT(2.0*(qbeame_exp(k)+qbeami_exp(k))/
     >     (sbeam_exp(k)*amassgas_exp*1.6726D-27))
c         mparsour_p=flow_beam_exp(k)*cosbeam*vbeam
c         mparsour_mm=flow_beam_exp(k-1)*cosbeam*vbeam
         if(itorque.eq.0)then
c torque denty in NT/M^2
           torque_exp(k)=(amassgas_exp*rmajor_exp*1.6726D-27)*
     >     cosbeam*vbeam*sbeam_exp(k)
         endif
          Mphi(k)= torque_exp(k)
          stress_tor_exp(k) = stress_tor_exp(k-1) +
     >     vprime(k,1)*torque_exp(k)*dr(k,1)
          dvoldr = vprime(k,2)*drhodr(k)*drhodr(k)
          tim = (ti_exp(k+1)+ti_exp(k))/2.0
          nem = (ne_exp(k+1)+ne_exp(k))/2.0
          fim = (fi_m(k+1)+fi_m(k))/2.0
          nim = fim*nem
          gradvexbm = (vexb_exp(k+1)-vexb_exp(k))/dr(k,2)
          ctorm = (c_tor(k+1)+c_tor(k))/2.0
          grad_c_tor=(c_tor(k+1)-c_tor(k))/dr(k,2)
          gradvphim =(ctorm*gradvexbm+grad_c_tor*vexbm)*cv
          eta_tor_exp(k) = -stress_tor_exp(k)*gradvphim/
     >    (1.62726D-8*amassgas_exp*nim*rmajor_exp*dvoldr
     >    *MAX(1.0D-12,gradvphim**2))
c         endif
         Mpar(k)=Mphi(k)*c_per(k)/c_tor(k)
c         write(*,'(i2,4x,1p2e12.4," Mphi")') k,rho(k),Mphi(k)
          stress_par_exp(k) = stress_tor_exp(k)*c_per(k)/c_tor(k)
      enddo
      eta_tor_exp(0)=eta_tor_exp(1)
      eta_tor_exp(mxgrid)=eta_tor_exp(mxgrid-1)
       Peaux_tot=volint(Peaux)
       Piaux_tot=volint(Piaux)
       Pohm_tot=volint(Pohpro)
       if(i_proc.eq.0)write(6,*)"Pohm = ",Pohm_tot       
       p_glob_exp = Peaux_tot+Piaux_tot+Pohm_tot
       if(i_proc.eq.0)write(6,*)"Paux = ",Peaux_tot,Piaux_tot
     >  ,Piaux_tot+Peaux_tot
       Peaux_tot=volint(Pe_alpha)
       Piaux_tot=volint(Pi_alpha)
       if(i_proc.eq.0)write(6,*)"P_alpha = ",Peaux_tot,Piaux_tot
     >  ,Piaux_tot+Peaux_tot
       if(i_proc.eq.0)write(6,*)"cmodel = ",cmodel
       if(i_proc.eq.0)write(6,*)"xparam_gf(10) = ",xparam_gf(10)
c
crew6.30
      do j=0,jmaxm
       pow_ei_cor_m(j)=0.
      enddo
c
c
c
c compute ExB shear on half grid using computed v_exb
c
      do k=1,mxgrid-1
         rminm=(rmin_exp(k+1)+rmin_exp(k))/2.0
         rhom=arho_exp*(rho(k+1)+rho(k))/2.0
         csdam=9.79D5*DSQRT(1.D3*(te_exp(k+1)+te_exp(k))/2.0)/
     >    (arho_exp*100.D0)/DSQRT(amassgas_exp)
         egamma_exp(k) = -cv/csdam*rminm/rhom*drhodr(k)*theta_exp(k)*
     >  (vexb_exp(k+1)-vexb_exp(k))/dr(k,2)
c       write(*,'i2,2x,0p6f10.5') k, rho(k), vexb_exp(k),
c     >       egamma_exp(k)
      enddo
c
c compute ExB shear on half grid using exp. Er
c Note: GYRO uses local Te in c_s/a and a=r(a)
c in GYRO while a=rho(a) in XPTOR
c Note #2: turn off smoothing of rmin,rmaj using average7_1d
c in datmap.f
c
      if(iexb.eq.2) then
c
      do k=1,mxgrid-1
        rminm=(rmin_exp(k+1)+rmin_exp(k))/2.D0
        rhom=arho_exp*(rho(k+1)+rho(k))/2.D0
c        if(igeo_tg.eq.0)then
c          b_unit = bt_exp
c        else
          b_unit = bt_exp*(rhom/rminm)*drhodr(k)
c        endif
        gradr_exp(k)=1.D0 / ( 1.D0 +drhodr(k)*
     >             ( rmaj_exp(k+1)-rmaj_exp(k) )/dr(k,2) )
        vexb_exp(k) = er_exp(k)/b_unit
c        write(*,'i2,2x,0p6f10.5') k, rho(k), er_exp(k),
c     >        b_unit, vexb_exp(k)
      enddo
c
c      write(*,*) 'arho = ',arho_exp, rmin_exp(mxgrid)
      do k=1,mxgrid-1
         rminm=(rmin_exp(k+1)+rmin_exp(k))/2.D0
         rhom=arho_exp*(rho(k+1)+rho(k))/2.D0
         qm=(q_exp(k+1)+q_exp(k))/2.D0
         temm=(te_exp(k+1)+te_exp(k))/2.D0
         gradrm=(gradr_exp(k+1)+gradr_exp(k))/2.D0
         csdam=9.79D5*DSQRT(1.D3*temm)/
     >    (arho_exp*100.D0)/DSQRT(amassgas_exp)
c         egamma_exp(k) = 1000.D0*rminm/rhom*drhodr(k)*
c     >     ((vexb_exp(k+1)-vexb_exp(k))/dr(k,2)
c     >     + zptheta_exp(k)*(vexb_exp(k+1)+vexb_exp(k))/2.D0)/
c     >     gradrm/csdam*
c     >     dsqrt(temm)/dsqrt(1.24416)
         egamma_exp(k) = 1000.D0*rminm/qm*drhodr(k)*
     >     ( (q_exp(k+1)/rmin_exp(k+1))*vexb_exp(k+1) -
     >       (q_exp(k)/rmin_exp(k))*vexb_exp(k) ) / dr(k,2) /
     >     gradrm/csdam*rmin_exp(mxgrid)/arho_exp*
     >     dsqrt(temm)/dsqrt(1.24416)
c       write(*,'i2,2x,0p6f10.5') k,rho(k),rmin_exp(k),vexb_exp(k),
c     >       vexb_m(k), gradr_exp(k), egamma_exp(k)
      enddo
      endif
c
      egamma_exp(0)=0.D0
      egamma_exp(mxgrid)=egamma_exp(mxgrid-1)
c
      j=0
c
      if(itport_pt(1).ne.0)then
        j=j+1
        do k=1,ngrid
          Told(j,k) = ne_m(k)
        enddo
      endif
      if(itport_pt(2).ne.0)then
        j=j+1
        do k=1,ngrid
           Told(j,k) = te_m(k)
        enddo
      endif
      if(itport_pt(3).ne.0)then
        j=j+1
        do k=1,ngrid
           Told(j,k) = ti_m(k)
        enddo
      endif
      if(itport_pt(4).ne.0)then
        j=j+1
        do k=1,ngrid
          Told(j,k) = ne_m(k)*vexb_m(k)
        enddo
      endif
      if(itport_pt(5).ne.0)then
        j=j+1
        do k=1,ngrid
          Told(j,k) = ne_m(k)*vpol_m(k)
        enddo
      endif
c
 150  format(2i4)
 310  format(0pe16.7,2x,0p301e20.11)
c
      return
      end
c
