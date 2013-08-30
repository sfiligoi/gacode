cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine glf2d_dv
c
      Implicit None
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      real*8 gradvphim
c
      call neoclassical
c
      if(ipert_gf.eq.0)then
        flowe_neo(jm) = vprime(jm,2)*neflux_neo
        flowi_neo(jm) = vprime(jm,2)*niflux_neo
        flowz_neo(jm) = vprime(jm,2)*nzflux_neo
        powe_neo(jm) = vprime(jm,2)*teflux_neo
        powi_neo(jm) = vprime(jm,2)*tiflux_neo
        powz_neo(jm) = vprime(jm,2)*tzflux_neo
        stress_tor_i_neo(jm) = vprime(jm,2)*vphiflux_neo
        stress_tor_z_neo(jm) = vprime(jm,2)*vphizflux_neo
        stress_par_i_neo(jm) = vprime(jm,2)*vparflux_neo
        stress_par_z_neo(jm) = vprime(jm,2)*vparzflux_neo
      endif
c
      if(imodel.eq.81)call glf23_dv
c
      if(imodel.eq.82)call tglf_dv
c
      nefluxm = neflux_neo + neflux_glf
      nifluxm = niflux_neo + niflux_glf
      nzfluxm = nzflux_neo + nzflux_glf
      tefluxm = teflux_neo + teflux_glf
      tifluxm = tiflux_neo + tiflux_glf
      tzfluxm = tzflux_neo + tzflux_glf
      vphifluxm = vphiflux_neo + vphiflux_glf
      vphizfluxm = vphizflux_neo + vphizflux_glf
      vparfluxm = vparflux_neo + vparflux_glf
      vparzfluxm = vparzflux_neo + vparzflux_glf
c
c save glf flows
c
      if(ipert_gf.eq.0)then
c save fluxes in MKS units
        flowe_glf(jm) = vprime(jm,2)*neflux_glf           !KA
        flowi_glf(jm) = vprime(jm,2)*niflux_glf           !KA
        flowz_glf(jm) = vprime(jm,2)*nzflux_glf           !KA
        powe_glf(jm) = vprime(jm,2)*teflux_glf            !MW
        powi_glf(jm) = vprime(jm,2)*tiflux_glf            !MW
        powz_glf(jm) = vprime(jm,2)*tzflux_glf            !MW
        stress_tor_i_glf(jm) = vprime(jm,2)*vphiflux_glf  !NT-M
        stress_tor_z_glf(jm) = vprime(jm,2)*vphizflux_glf !NT-M
        stress_par_i_glf(jm) = vprime(jm,2)*vparflux_glf  !NT-M
        stress_par_z_glf(jm) = vprime(jm,2)*vparzflux_glf !NT-M
      endif
c
      if(imodel.eq.81.or.imodel.eq.82)then
        call adhoc_dv
c
        if(ipert_gf.eq.0)then
        flowe_adhoc(jm) = vprime(jm,2)*neflux_adhoc
        flowi_adhoc(jm) = vprime(jm,2)*niflux_adhoc
        flowz_adhoc(jm) = vprime(jm,2)*nzflux_adhoc
        powe_adhoc(jm) = vprime(jm,2)*teflux_adhoc
        powi_adhoc(jm) = vprime(jm,2)*tiflux_adhoc
        powz_adhoc(jm) = vprime(jm,2)*tzflux_adhoc
        stress_tor_i_adhoc(jm) = vprime(jm,2)*vphiflux_adhoc
        stress_tor_z_adhoc(jm) = vprime(jm,2)*vphizflux_adhoc
        stress_par_i_adhoc(jm) = vprime(jm,2)*vparflux_adhoc
        stress_par_z_adhoc(jm) = vprime(jm,2)*vparzflux_adhoc
        endif

c
        nefluxm = nefluxm + neflux_adhoc
        nifluxm = nifluxm + niflux_adhoc
        nzfluxm = nzfluxm + nzflux_adhoc
        tefluxm = tefluxm + teflux_adhoc
        tifluxm = tifluxm + tiflux_adhoc
        tzfluxm = tzfluxm + tzflux_adhoc
        vphifluxm = vphifluxm + vphiflux_adhoc
        vphizfluxm = vphizfluxm + vphizflux_adhoc
        vparfluxm = vparfluxm + vparflux_adhoc
        vparzfluxm = vparzfluxm + vparzflux_adhoc
      endif
c
c      if(imodel.eq.2) call ifsppl_dv
c
      if(imodel.eq.7) call mixshear_dv
c
c      if(imodel.eq.99) call chi1_dv
c
      if(imodel.eq.-1) call test_dv
c
      if(ipert_gf.eq.0)then
c save fluxes in MKS units
        flowe_m(jm) = vprime(jm,2)*nefluxm           !KA
        flowi_m(jm) = vprime(jm,2)*nifluxm           !KA
        flowz_m(jm) = vprime(jm,2)*nzfluxm           !KA
        powe_m(jm) = vprime(jm,2)*tefluxm            !MW
        powi_m(jm) = vprime(jm,2)*tifluxm            !MW
        powz_m(jm) = vprime(jm,2)*tzfluxm            !MW
        stress_tor_i_m(jm) = vprime(jm,2)*vphifluxm  !NT-M
        stress_tor_z_m(jm) = vprime(jm,2)*vphizfluxm !NT-M
        stress_par_i_m(jm) = vprime(jm,2)*vparfluxm  !NT-M
        stress_par_z_m(jm) = vprime(jm,2)*vparzfluxm !NT-M
      endif
c
c      summ the ions
c
      tifluxm = tifluxm + tzfluxm
      vphifluxm = vphifluxm + vphizfluxm
      vparfluxm = vparfluxm + vparzfluxm
c
      if(ipert_gf.eq.0)then
       diffgb_m(jm) = -gradnem*nefluxm
     >  /(cgyrobohm_m(jm)*1.6022D-3*MAX(1.0D-10,gradnem*gradnem))
       chiegb_m(jm) = -gradtem*tefluxm
     >  /(cgyrobohm_m(jm)*1.6022D-3*nem*MAX(1.0D-10,gradtem*gradtem))
       chiigb_m(jm) = -gradtim*tifluxm
     >  /(cgyrobohm_m(jm)*1.6022D-3*nim*MAX(1.0D-10,gradtim*gradtim))
       gradvphim = cv*(vphi_m(jm+1)-vphi_m(jm))/dr(jm,2)
       etagb_phi_m(jm) = -gradvphim*vphifluxm
     >  /(cgyrobohm_m(jm)*1.6726D-8*amassgas_exp*nim
     >     *MAX(1.0D-10,gradvphim*gradvphim))
       chiegb_etg_m(jm) = -gradtem*tefluxm_etg
     >  /(cgyrobohm_m(jm)*1.6022D-3*nem*MAX(1.0D-10,gradtem*gradtem))
       chiineogb_m(jm) = -gradtim*(tiflux_neo+tzflux_neo)
     >  /(cgyrobohm_m(jm)*1.6022D-3*nim*MAX(1.0D-10,gradtim*gradtim)) 
       chieneogb_m(jm) = -gradtem*teflux_neo
     >  /(cgyrobohm_m(jm)*1.6022D-3*nem*MAX(1.0D-10,gradtem*gradtem)) 
      else
c this is needed in order for mpi to sum up right
        diffgb_m(jm) = 0.0
        chiegb_m(jm) = 0.0
        chiigb_m(jm) = 0.0
        etagb_phi_m(jm) = 0.0
        chiegb_etg_m(jm) = 0.0
        chiineogb_m(jm) = 0.0
        chieneogb_m(jm) = 0.0
        rhosda_m(jm) = 0.0
        csda_m(jm) = 0.0
        cgyrobohm_m(jm) = 0.0
      endif
c
      return
      END 
c
      SUBROUTINE chi1_dv
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      real*8 gfac
c
      cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*DSQRT(tem*1.D3)/(arho_exp*100.D0)
     >  *(1.02D2*DSQRT(tem*1.D3)/bt_exp/1.D4)**2*DSQRT(amassgas_exp)
      rhosda_m(jm)=((1.02D2*DSQRT(tem*1.D3))/bt_exp/1.D4)
     >  *DSQRT(amassgas_exp)/(arho_exp*100.D0)
      zpmni = -arho_exp*gradnim/nim
      zpmte=-arho_exp*gradtem/tem
      zpmti=-arho_exp*gradtim/tim
      zpni_m(jm)=zpmni
      zpte_m(jm)=zpmte
      zpti_m(jm)=zpmti
      csdam=9.79D5*DSQRT(tem*1.D3)/
     &      (arho_exp*100.D0)/DSQRT(amassgas_exp)
      csda_m(jm)=csdam
      zpmni=zpmne
c
       if(igeo_m.ge.1) gfac=geofac(jm)
c
       chiineo_m=0.D0
       chieneo_m=0.D0
       chiineogb_m(jm)=0.D0
       chieneogb_m(jm)=0.D0
c
       diffnem=0.D0
       chietem=1.D0
       chiitim = 1.D0
       chie_m(jm)=chietem
       chii_m(jm)=chiitim
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
c... compute fluxes, power flows
c
      tefluxm = tem*nem*zpmte*chietem/arho_exp
      tifluxm = tim*nim*zpmti*chiitim/arho_exp
      powem = kevdsecpmw*tem*nem*1.D19/arho_exp*gradrhosq_exp(jm)*
     &        sfactor(jm)*(chietem*zpmte)
     &        +xconv*1.5D0*tem*flow_exp(jm)
      powim = kevdsecpmw*tim*nim*1.D19/arho_exp*gradrhosq_exp(jm)*
     &        sfactor(jm)*(chiitim*zpmti)
     &        +xconv*1.5D0*tim*flow_exp(jm)
      powem = 1.6022D-3*vprime(jm,2)*tefluxm + pow_ei_cor_m(jm)
      powim = 1.6022D-3*vprime(jm,2)*tifluxm - pow_ei_cor_m(jm)
c
      return
      END 
c
      SUBROUTINE mixshear_dv
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      real*8 ceqei_ms, chiemi_ms, facbr_ms, facgb_ms, fshea_ms,
     &       dbohm_ms, gyro_ms, rlte_ms, rlti_ms, gfac
c
       cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*DSQRT(tem*1.D3)/(arho_exp*100.D0)
     >  *(1.02D2*DSQRT(tem*1.D3)/bt_exp/1.D4)**2*DSQRT(amassgas_exp)
       rhosda_m(jm)=((1.02D2*DSQRT(tem*1.D3))/bt_exp/1.D4)
     >  *DSQRT(amassgas_exp)/(arho_exp*100.D0)
       csdam=9.79D5*DSQRT(tem*1.D3)/
     &      (arho_exp*100.D0)/DSQRT(amassgas_exp)
       csda_m(jm)=csdam
c
       ceqei_ms  = 3.5D0
       chiemi_ms = 5.D-02
       facbr_ms  = 8.61D-03
       facgb_ms  = 5.07D-01
c
       if ( shat_exp(jm) .le. 0. ) then
         fshea_ms = 0.D0
       else
         fshea_ms = shat_exp(jm)**2.D0/(1.D0+shat_exp(jm)**3.D0)
       endif
c
       dbohm_ms = tem*1000.D0/(16.D0*bt_exp)
       gyro_ms  = 1.0D-4*dsqrt(tem*1000.D0)/bt_exp
       zpmni = -arho_exp*gradnim/nim
       zpmte=-arho_exp*gradtem/tem
       zpmti=-arho_exp*gradtim/tim
       rlte_ms  = - arho_exp*gradtem/tem
       rlti_ms  = - arho_exp*gradtim/tim
       gfac=1.D0
       if(igeo_m.ge.1) gfac=geofac(jm)

       chietem=cmodel_e*gfac
     &       *abs( facbr_ms*rlte_ms*dbohm_ms*q_exp(jm)**2.D0
     &       *fshea_ms + facgb_ms*rlte_ms*dbohm_ms*gyro_ms/arho_exp )
       chiitim=cmodel_i*gfac
     &       *abs( ceqei_ms*facbr_ms*rlti_ms*dbohm_ms*q_exp(jm)**2.D0
     &       *fshea_ms + facgb_ms*rlti_ms*dbohm_ms*gyro_ms/arho_exp )
       diffnem=0.D0
c
       chietem=dmax1(chietem,chiemi_ms)
c
       chie_m(jm)=chietem
       chii_m(jm)=chiitim
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
c... compute fluxes, power flows
c
      tefluxm = tem*nem*zpmte*chietem/arho_exp
      tifluxm = tim*nim*zpmti*chiitim/arho_exp
      powem = kevdsecpmw*tem*nem*1.D19/arho_exp*gradrhosq_exp(jm)*
     &        sfactor(jm)*(chietem*zpmte)
     &        +xconv*1.5D0*tem*flow_exp(jm)
      powim = kevdsecpmw*tim*nim*1.D19/arho_exp*gradrhosq_exp(jm)*
     &        sfactor(jm)*(chiitim*zpmti)
     &        +xconv*1.5D0*tim*flow_exp(jm)
      powem = 1.6022D-3*vprime(jm,2)*tefluxm + pow_ei_cor_m(jm)
      powim = 1.6022D-3*vprime(jm,2)*tifluxm - pow_ei_cor_m(jm)
c
      return
      END 
c
      SUBROUTINE glf23_dv
c
c... imodel=0 for GLF23 model
c      
c
      Implicit None
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      integer i,j,k,ik
      real*8 exchgb_local,gfac,diffgb_local
      real*8 fc,vnewk3x,alpha_neo_hold,akappa1
      real*8 aiwt_jp1,xnimp_jp1,xnimp
      real*8 rmajm,rminm,rhom
      real*8 fac_imp_flow
      real*8 zeffm,qm, mass_factor
      real*8 vneom(nspecies),vdiam(nspecies)
      real*8 gradvneom(nspecies),gradvdiam(nspecies)
      real*8 gradvphim
      real*8 cparm,cperm,ctorm,grad_c_par,grad_c_per,grad_c_tor
      real*8 apolm,atorm,grad_a_pol,grad_a_tor
      real*8 chietem_etg
c
      cxnu = 1.D0
      qm =(q_exp(jm+1)+q_exp(jm))/2.D0
      rmajm=(rmaj_exp(jm+1)+rmaj_exp(jm))/2.D0
      rhom=arho_exp*(rho(jm+1)+rho(jm))/2.D0
      rminm=(rmin_exp(jm+1)+rmin_exp(jm))/2.D0
      zeffm=(zeff_exp(jm+1)+zeff_exp(jm))/2.D0
      zpmne = -arho_exp*gradnem/nem
      if(zpmne.eq.0.D0)zpmne=1.0D-7
      zpmni = -arho_exp*gradnim/nim
      if(zpmni.eq.0.D0)zpmni=1.0D-7
      zpmte = -arho_exp*gradtem/tem
      if(zpmte.eq.0.D0)zpmte=1.0D-7
      zpmti = -arho_exp*gradtim/tim
      if(zpmti.eq.0.D0)zpmti=1.0D-7
      csdam=9.79D5*DSQRT(tem*1.D3)/
     > (arho_exp*100.D0)/DSQRT(amassgas_exp)
      csda_m(jm) = csdam
      vdiam(2) = (vdia_m(2,jm+1)+vdia_m(2,jm))/2.0
      vneom(2) = (vneo_m(2,jm+1)+vneo_m(2,jm))/2.0
      cparm = (c_par(jm+1)+c_par(jm))/2.0
      cperm = (c_per(jm+1)+c_per(jm))/2.0
      ctorm = (c_tor(jm+1)+c_tor(jm))/2.0
      apolm = (a_pol(jm+1)+a_pol(jm))/2.0
      atorm = (a_tor(jm+1)+a_tor(jm))/2.0
      gradvdiam(2) =(vdia_m(2,jm+1)-vdia_m(2,jm))/dr(jm,2)
      gradvneom(2) =(vneo_m(2,jm+1)-vneo_m(2,jm))/dr(jm,2)
      grad_c_par=(c_par(jm+1)-c_par(jm))/dr(jm,2)
      grad_c_per=(c_per(jm+1)-c_per(jm))/dr(jm,2)
      grad_c_tor=(c_tor(jm+1)-c_tor(jm))/dr(jm,2)
      grad_a_pol=(a_pol(jm+1)-a_pol(jm))/dr(jm,2)
      grad_a_tor=(a_tor(jm+1)-a_tor(jm))/dr(jm,2)
c
c      gamma_p_m(jm) = -(cv/csdam)*drhodr(jm)*
c     >  (apolm*(gradvpolm+gradvneom(2))+atorm*(gradvexbm+gradvdiam(2))
c     >  +(vpolm+vneom(2))*grad_a_pol+(vexbm+vdiam(2))*grad_a_tor)
c      egamma_m(jm) = -cv/csdam*(rminm/rhom)*drhodr(jm)
c     >  *theta_exp(jm)*gradvexbm
c      if(jm.eq.ngrid-1.and.itport_pt(5).eq.0)egamma_m(jm)=0.0
c       gradvphim = -gamma_p_m(jm)*csdam/cv
c       gradvphim =
c     >  cperm*(gradvpolm+gradvneom(2))+ctorm*(gradvexbm+gradvdiam(2))
c     > +grad_c_per*(vpolm+vneom(2))+grad_c_tor*(vexbm+vdiam(2))
c
        gamma_e_gf = -cv/csdam*(rminm/rhom)*drhodr(jm)
     >  *theta_exp(jm)*gradvexbm
      if(ipert_gf.eq.0)egamma_m(jm)=csdam*gamma_e_gf/cv
      if(iexb.ge.1)then
c        gamma_e_gf = egamma_exp(jm)
        gamma_e_gf =(cv/csdam)* egamma_exp(jm)
      endif
      if(doppler_shear_model.eq.1)then
c        gamma_e_gf = -cv/csdam*(rminm/rhom)*drhodr(jm)
c     >  *theta_exp(jm)*(vexb_m(jm+1)-vexb_m(jm))/dr(jm,2)
        gamma_e_gf = (cv/csdam)*doppler_shear_m(jm)
      endif
      gamma_p_gf = -(cv/csdam)*drhodr(jm)*
     > (gradvexbm+gradvdiam(2)+(apolm/atorm)*(gradvpolm+gradvneom(2)))
c     >  +(vpolm+vneom(2))*grad_a_pol+(vexbm+vdiam(2))*grad_a_tor)
      if(jm.eq.ngrid-1.and.itport_pt(5).eq.0)gamma_e_gf=0.0
      if(ipert_gf.eq.0)gamma_p_m(jm) = gamma_p_gf
      exch_gf=0.D0
c   local rho_star
      if(bt_flag.gt.0)then
c    use effective B-field
       rhosda_m(jm)=((1.02D2*DSQRT(tem*1.D3))/bteff_exp(jm)
     > /1.D4)*DSQRT(amassgas_exp)/(arho_exp*100.D0)
      else
       rhosda_m(jm)=((1.02D2*DSQRT(tem*1.D3))/bt_exp/1.D4)
     >  *DSQRT(amassgas_exp)/(arho_exp*100.D0)
      endif
c   local gyrobohm unit of diffusion
      if(bt_flag.ge.2)then
c    use effective B-field
       cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*DSQRT(tem*1.D3)/(arho_exp*100.D0)
     >  *(1.02D2*DSQRT(tem*1.D3)/bteff_exp(jm)/1.D4)**2
     >  *DSQRT(amassgas_exp)
      else
c   local gyrobohm unit of diffusion
       cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*DSQRT(tem*1.D3)/(arho_exp*100.D0)
     >  *(1.02D2*DSQRT(tem*1.D3)/bt_exp/1.D4)**2*DSQRT(amassgas_exp)
      endif 
       betae_m(jm) = 4.03D-3*nem*tem/(bt_exp**2)
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
       vnewk3x=
     >   0.117D0*nem/DSQRT(tem)**3/DSQRT(tim)*(arho_exp)*
     >   DSQRT(amassgas_exp/2.D0)
c
       xnu_m(jm) =vnewk3x*DSQRT(2.D0*tim/tem)
c
       if(ipert_gf.eq.0)alpha_m(jm)=drhodr(jm)*qm**2*rmajm
     >   *betae_m(jm)*((tim*nim/tem/nem)*
     >   (zpmni+zpmti)+zpmne+zpmte)/arho_exp
c bananna regime ie collisionless limit formulas
        alpha_neo_hold=alpha_neo
        fc=1-1.46D0*DSQRT(rminm/rmajm)+
     >      0.46D0*DSQRT(rminm/rmajm)**3
        akappa1=0.8839D0*fc/(0.3477D0+0.4058D0*fc)
        if(irot1.eq.1) alpha_neo=-akappa1+1.D0

       rmaj_gf=rmajm/arho_exp
       rmin_gf=rminm/arho_exp
       q_gf=qm
       betae_gf=cbetae*betae_m(jm)
       shat_gf=shat_exp(jm)
        if(igeo_m.ge.2) shat_gf=shat_exp(jm)*drhodrrrho(jm)
       alpha_gf=xalpha*alpha_exp(jm)
       if(ialphastab.gt.0) then
        alpha_gf=xalpha*alpha_m(jm)
       endif
       elong_gf=elong_exp(jm)
       if(igeo_m.eq.-1) elong_gf=elonga_exp
       xnu_gf=cxnu*xnu_m(jm)
c       write(*,*) jm, rho(jm), xnu_gf, 'xnu'
       taui_gf=tim/tem
       amassgas_gf=amassgas_exp
       rlte_gf=zpmte
       rlti_gf=zpmti
       rlne_gf=zpmne
       rlni_gf=rlne_gf
       if(i_dengrad.ge.1)then
         rlni_gf=zpmni
       endif
       dil_gf=0.D0
       apwt_gf=1.D0
       aiwt_gf=0.D0
       rlnimp_gf=1.D0
       zpmnimp=1.D0
       if(i_dengrad.eq.2) dil_gf=1.D0-nim/nem
       if(i_dengrad.eq.3) then
cgms not hooked up right
        apwt_gf=nim/nem
        aiwt_jp1=(zeffm*nem-ni_exp(jm+1)
     >    -nfast_exp(jm+1))/(zimp_gf**2*ne_exp(jm+1))
        xnimp_jp1=aiwt_jp1*ne_exp(jm+1)
        aiwt_gf=(zeffm*nem-nim
     >    -nfast_exp(jm))/(zimp_gf**2*ne_exp(jm))
        xnimp=aiwt_gf*ne_exp(jm)
        
        zpmnimp=-(xnimp_jp1-xnimp)/(rho(jm+1)-rho(jm))/
     >  ((xnimp_jp1+xnimp)/2.0D0)
        rlnimp_gf=zpmnimp
c zimp_gf and amassimp_gf are inputs
       endif
       if(igeo_m.ne.0)then
        rlte_gf=zpmte*DSQRT(elong_exp(jm))
        rlti_gf=zpmti*DSQRT(elong_exp(jm))
        rlne_gf=zpmne*DSQRT(elong_exp(jm))
        rlni_gf=rlne_gf
        if(i_dengrad.ge.1) then
         rlne_gf=zpmne*DSQRT(elong_exp(jm))
         rlni_gf=zpmni*DSQRT(elong_exp(jm)) 
        endif
        rlnimp_gf=zpmnimp*DSQRT(elong_exp(jm))
       endif
        if(igeo_m.ge.1) then
         rlte_gf=zpmte*drhodr(jm)
         rlti_gf=zpmti*drhodr(jm)
         rlne_gf=zpmne*drhodr(jm)
         rlni_gf=rlne_gf
         rlnimp_gf=zpmnimp*drhodr(jm)
c
         if(i_dengrad.ge.1) then
          rlne_gf=zpmne*drhodr(jm)
          rlni_gf=zpmni*drhodr(jm)
         endif
        endif
        if(igeo_m.eq.-1) then
         rlte_gf=zpmte*DSQRT(elonga_exp)
         rlti_gf=zpmti*DSQRT(elonga_exp)
         rlnimp_gf=zpmnimp*DSQRT(elonga_exp)
         rlne_gf=zpmne*DSQRT(elonga_exp)
         rlni_gf=rlne_gf
         if(i_dengrad.ge.1) then
          rlne_gf=zpmne*DSQRT(elonga_exp)
          rlni_gf=zpmni*DSQRT(elonga_exp)
          rlnimp_gf=zpmnimp*DSQRT(elonga_exp)
         endif
       endif
ctemp
c       rlni_gf=rlne_gf
ctemp
       gamma_star_gf=0.0D0
       gamma_mode_gf=0.0D0
c
 11       format(2x,1pe12.5)
c
c       write(*,*) 'q,shat = ',q_gf,shat_gf
c       write(*,*) 'a/Lt = ',rlti_gf,rlte_gf
c       write(*,*) 'a/Ln = ',rlni_gf,rlne_gf
c       write(*,*) 'R/a,r/a = ',rmaj_gf,rmin_gf
       call glf2d
c switch from ion particle diffusivity to electron particle diffusivity
       diff_gf = nim*rlni_gf*diff_gf/(nem*rlne_gf)
c
c       write(*,*) 'gamma_gf = ',jm,gamma_gf(1)
c
c 4/25/96 note: chie_gf and chii_gf from glf2d are energy diffusivies
c but here we can redefine them as heat diffusivites if we
c take the convection from experiment using xconv=1 or 5./3.
cgms disable xconv
c       chie_gf=chie_gf-xconv*1.5D0*rlni_gf/rlte_gf*diff_gf
c       chii_gf=chii_gf-xconv*1.5D0*rlni_gf/rlti_gf*diff_gf
       ky_j(jm)=xkyf_gf
       gamma_j(jm,1)=gamma_gf(1)
       gamma_j(jm,2)=gamma_gf(2)
       gamma_j(jm,3)=gamma_gf(3)
       gamma_j(jm,4)=gamma_gf(4)
c 
       freq_j(jm,1)=freq_gf(1)
       freq_j(jm,2)=freq_gf(2)
       freq_j(jm,3)=freq_gf(3)
       freq_j(jm,4)=freq_gf(4)
c
       phi_norm_j(jm,1)=phi_norm_gf(1)
       phi_norm_j(jm,2)=phi_norm_gf(2)
       phi_norm_j(jm,3)=phi_norm_gf(3)
       phi_norm_j(jm,4)=phi_norm_gf(4)
c 
       do ik=1,ikymax_gf
        gamma_k_j(ik,jm)=gamma_k_gf(1,ik)
        freq_k_j(ik,jm) =freq_k_gf(1,ik)
        chie_k_j(ik,jm) = chie_k_gf(ik)
        chii_k_j(ik,jm) = chii_k_gf(ik)
       enddo
c
       if(ipert_gf.eq.0)then
         anrate_m(jm)=gamma_gf(1)
c         write(*,*)jm,anrate_m(jm)
         dnrate_m(jm)=0.D0
         dtnrate_m(jm)=0.D0
         anfreq_m(jm)=freq_gf(1)
         dnfreq_m(jm)=0.D0
       endif
c
        gfac=1.D0
        if(igeo_m.ge.1) gfac=geofac(jm)
        if(igeo_m.eq.-1) 
     >   gfac=(1.D0+elonga_exp**2)/(2.D0*elonga_exp)/gradrhosq_exp(jm)
cgms  rescale gfac to agree with shooting code
cgms       gfac=gfac*gradrhosq_exp(jm)*sfactor(jm)/vprime(jm,2)
       chiegb_m(jm)=chie_gf*cmodel
       chiigb_m(jm)=chii_gf*cmodel
       diffgb_m(jm)=diff_gf*cmodel
       diffgb_im_m(jm)=diff_im_gf*cmodel
c note diff_im_gf used only for diagnostics
c  to be compared with diff_gf or chii_gf
c
       exchgb_m(jm)=exch_gf*cmodel
       etagb_phi_m(jm)=eta_phi_gf*cmodel
       etagb_par_m(jm)=eta_par_gf*cmodel
       etagb_per_m(jm)=eta_per_gf*cmodel   
       chiegb_etg_m(jm)=chie_e_gf*cmodel
c  
c exch_m in MW/m**3  
c
cgms      exch_m(jm)=1.D19*
cgms     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel 
c
c exch_m is directly related to the flow
c for a single mode branch exch_gf=-(-freq_gf(1)/xkyf_gf)*diff_gf*rln_gf.
c we can  not expect to find exch_m without knowing flow_exp as input.
c and solving self consistant flow eq. flown=flow_exp for density
c density profile.
c
c however, knowing freq_gf(1) from the gf model we can compute exch_exp  
c from flow_exp using
c       flowm=kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*
c     >       sfactor(jm)*(difftem*zpmte+difftim*zpmti+diffnem*zpmne)
c we have:

cgms       diffgb_local=flow_exp(jm)/
cgms     > (kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*sfactor(jm)*
cgms     > zpmne_q)/cgyrobohm_m(jm)
c
cgms       exchgb_local=-(-freq_gf(1)/xkyf_gf)*diffgb_local*rlni_gf
c 
cgms       exch_exp(jm)=1.e19*
cgms     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exchgb_local
c   
c       exch_exp(jm)=flow_exp(jm)*tem*(-1.)*(-freq_gf(1)/xkyf_gf)*
c     > DSQRT(elong_exp(jm))*arho_exp/gradrhosq_exp(jm)/sfactor(jm)
c     >/arho_exp(jm)**2 
c
c   note electron(ion) wave freq > 0(<0) cool(heat) electrons
c (-1) denotes electron to ion
c
c to emphasize, we can not know exch_exp better than we know flow_exp
c
c chietem, chiitim, diffen in m**2/sec
c
      chietem=cmodel*gfac*chie_gf*cgyrobohm_m(jm)
      chietem=chietem*cmodel_e
cgms      chietim=0.
cgms       chienem=0.
cgms       chiitem=0.
      chiitim=cmodel*gfac*chii_gf*cgyrobohm_m(jm)
      chiitim=cmodel_i*chiitim
cgms       chiinem=0. 
cgms       difftem=0.
cgms       difftim=0.
c
      diffnem=cmodel*gfac*diff_gf*cgyrobohm_m(jm)      
      etaparm=cmodel*gfac*eta_par_gf*cgyrobohm_m(jm)
      etaexbm=cmodel*gfac*eta_per_gf*cgyrobohm_m(jm)
      etaphim=cmodel*gfac*eta_phi_gf*cgyrobohm_m(jm)
      etaphim=DABS(etaphim)
      etaparm=DABS(etaparm)
      chietem_etg = cmodel*gfac*chie_e_gf*cgyrobohm_m(jm)
c      etaexbm=etaparm
cgms temporary
c       etaexbm = DABS(chiitim)
c       etaparm = etaexbm
c save coeficients
      diff_m(jm)=diffnem
      chie_m(jm)=chietem
      chii_m(jm)=chiitim
      etapar_m(jm)=etaparm
      etaexb_m(jm)=etaexbm
      etaphi_m(jm)=etaphim
c
      chiegb_m(jm)=chietem/cgyrobohm_m(jm)
      chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
      diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
      etagb_phi_m(jm)=etaphim/cgyrobohm_m(jm)
      etagb_par_m(jm)=etaparm/cgyrobohm_m(jm)
      etagb_per_m(jm)=etaperm/cgyrobohm_m(jm)
      etagb_phi_m(jm)=etaphim/cgyrobohm_m(jm)
      chiegb_etg_m(jm)=chietem_etg/cgyrobohm_m(jm)
c
c compute fluxes
c
c switch to electron particle flux
       neflux_glf = (1.6022D-3)*nem*zpmne*diffnem/arho_exp
       teflux_glf = (1.6022D-3)*tem*nem*zpmte*chietem/arho_exp
       tiflux_glf = (1.6022D-3)*tim*nim*zpmti*chiitim/arho_exp
c use diagonal form for agreement with previous versions
       vphiflux_glf =(1.6726D-8) 
     >   *rmajor_exp*amassgas_exp*nim*etaphim*gamma_p_gf*csdam
       vparflux_glf = vphiflux_glf*c_per(jm)/c_tor(jm)
c
       mass_factor = 1.0 +amassimp_exp*fz_m(jm)/(amassgas_exp*fi_m(jm))
       vphiflux_glf = vphiflux_glf*mass_factor
       vparflux_glf = vparflux_glf*mass_factor
       niflux_glf = neflux_glf
       nzflux_glf = 0.0
       tzflux_glf = 0.0
       vphizflux_glf = 0.0
       vparzflux_glf = 0.0
       tefluxm_etg = (1.6022D-3)*tem*nem*zpmte*chietem_etg/arho_exp
c      
 50   format(2x,i2,2x,0p1f10.6,1p6e12.4)
c
      return
      END   !SUBROUTINE glf23_dv
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE neoclassical
c
      include '../inc/model.m'
c
      if(use_xneo_m.le.2)then
        CALL get_kapisn
      else
        CALL get_neo
      endif
c
      RETURN
      END SUBROUTINE neoclassical
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      SUBROUTINE get_kapisn
c
c... Calculate neoclassical chii from modelled profiles
c... Use adapted version KAPISN_ONE that calculats chii at one position
c... instead of the full profile
c... Note that arrays are name.nc while single parameters are name.nco
c
c
      use tglf_tg
      Implicit None
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      real*8 ztau_nco, zrho_nco, rhosdam, rhom, rminm,
     >       rmajm,b_unit
      real*8 taui,rhoi,lnlamda
      real*8 wthi,nu_star,d_zeff,zeffm,g_ft
      real*8 fcm,r_eps,nu_plateau,nu_norm
      real*8 nu_braginski,nu_banana,vthi,x13
      real*8 kpol1,kpol2,cnc
      real*8 thetam,qm
      real*8 ctorm, grad_c_tor
c
c      a_unit_exp = rmin_exp(mxgrid)
c      if(igeo_tg.eq.0)a_unit_exp=arho_exp
      rhom=arho_exp*(rho(jm+1)+rho(jm))/2.D0
      rminm=(rmin_exp(jm+1)+rmin_exp(jm))/2.D0
      rmajm=(rmaj_exp(jm+1)+rmaj_exp(jm))/2.0
      drhodr(jm) = arho_exp*(rho(jm+1)-rho(jm))/
     >             (rmin_exp(jm+1)-rmin_exp(jm))
      if(igeo_m.eq.0)then
        b_unit = bt_exp
c        if(bt_flag.gt.0)b_unit = bteff_exp(jm)
      else
        b_unit = bt_exp*(rhom/rminm)*drhodr(jm)
      endif
      csdam=9.79D5*DSQRT(tem*1.D3)/
     > (a_unit_exp*100.D0)/DSQRT(amassgas_exp)
      rhosdam=(1.02D2*DSQRT((tem*1.D3)*amassgas_exp)/
     >  (b_unit*1.D4))/(a_unit_exp*100.D0)
      cgyrobohm_m(jm)= csdam*(rhosdam*a_unit_exp)**2
c
c set all chanels to zero
c
      deneo_m(jm) = 0.0
      chieneo_m(jm) = 0.0
      chiineo_m(jm) = 0.0
      etaphineo_m(jm) = 0.0
c
      if (use_xneo_m.eq.0) then
        chieneo_m(jm)=0.D0 ! no electron chi from kapisn module
        chiineo_m(jm)=chiineo_exp(jm)
        chiineogb_m(jm)=chiineogb_exp(jm)
      elseif (use_xneo_m.eq.1) then
        ng_nc=1                     !number of hyd. species
        aplasm_nc(1)=amassgas_exp   !array of atomic masses of hyd. species
c        btf_nc=bt_exp               !tf (tesla) fld at cntr of outer flux
        btf_nc = b_unit
        drshaf_nc=0.0D0             !shaf shift of outer zone bdr. (cm)
cgms        numzones_nc=nj_d           !number of zones
cgms        rminor_nc=0.5D0*(rmin_exp(jm)+rmin_exp(jm+1))*1.D2  !on half grid
cgms        rmajor_nc=rmaj_exp(mxgrid)*1.D2
        rminor_nc = rminm*1.D2
        rmajor_nc = rmajm*1.D2
c
          rhoel_nco=dabs(nem*1.D13)             !electron density (cm**-3)
          rhob_nco(1)=dabs(nem*fi_m(jm)*1.D13)  !hyd. spec. den s (cm**-3)
          rhi_nco=dabs(nem*fz_m(jm)*1.D13)      !z.c. of av. impurity s density  
          rhoi_nco=rhi_nco+rhob_nco(1)   !z.c. of total ion density
          zeff_nco=(zeff_exp(jm+1)+zeff_exp(jm))/2.D0 !z.c. of plasma zeff
        if(nkimod_nc.eq.4)then  ! pure plasma CH model with b_unit 
          rhob_nco(1)=rhoel_nco
          rhi_nco = 0.0
          rhoi_nco = rhoel_nco
          zeff_nco = 1.0
          btf_nc = b_unit
          drshaf_nc = (rmaj_exp(jm+1)-rmaj_exp(jm))
     >              /(rmin_exp(jm+1)-rmin_exp(jm))
        endif
          te_nco=dabs(tem*1.D3)          !z.c.of Te (ev)
          ti_nco=dabs(tim*1.D3)          !z.c. of Ti (ev)
          q_nco=(q_exp(jm+1)+q_exp(jm))/2.D0          !z.c. of safety factor
          aimp_nco=pimpa
          xzimp_nco=pimpz
c
c       if(ipert_gf.eq.0.and.jm.eq.10)then
c        write(6,*)rmin_exp(mxgrid),rminavnpsi_d(nj_d)
c        write(6,*)"KAPISN input from glf2d_dv"
c        write(6,*)ng_nc,numzones_nc,btf_nc,drshaf_nc,rminor_nc,rmajor_nc
c        write(6,*)ng_nc,aplasm_nc(1),rhoel_nco,rhob_nco
c        write(6,*)rhoi_nco,te_nco,ti_nco,zeff_nco,q_nco
c        write(6,*)aimp_nco,xzimp_nco,rhi_nco
c       endif
c
c...    * means: default value given in input.mlt/mlt0in
c...    Note that the comments might be screwed up since some commentlines
c...    had to be deleted otherwise BASIS could not handle it (...) 
c...    See file neoclass.basf for the unscrewed comments
c
        call kapisn_one(
     >          nkimod_nc,         !*kapai model nr desired
     >          aimp_nco,          !*atomic mass of av. impurity
     >          xzimp_nco,         !*atomic number of av. impurity
     >          aplasm_nc,         !array of atomic masses of hyd. species
     >          ng_nc,             !number of hyd. species
     >          rhoel_nco,         !zone centered electron density (cm**-3)
     >          rhob_nco,          !z.c. array of hyd. spec. den s (cm**-3)
     >          rhi_nco,           !z.c. array of av. impurity s density  
     >          rhoi_nco,          !z.c. array of total ion density       
     >          te_nco,            !z.c. array of Te (ev)
     >          ti_nco,            !z.c. array of Ti (ev)
     >          zeff_nco,          !z.c. array of plasma zeff
     >          q_nco,             !z.c. array of safety factor
     >          btf_nc,            !tf (tesla) fld at cntr of outer flux 
     >          drshaf_nc,         !shaf shift of outer zone bdr. (cm)
     >          rminor_nc,         !plasma minor radius local on half grid (cm)
     >          rmajor_nc,         !major radius (cntr of outer flux) (cm)
     >          istringer_nc,      !Stringer correction
     >          xkapi_nco,         !o Neo-Class Ion thermal diff. (cm**2/sec)
     >          xnstari_nco,       !o nu-star-ions, 
     >          xnstare_nco,       !o nu-star-elecs, 
     >          ztau_nco,          !o ion collision time
     >          zrho_nco,          !o ion poloidal gyro-radius
     >          zfluxlim_nco)      !o flux lim flow max temp grad length
c
c      if(ipert_gf.eq.0.and.jm.eq.10)then
c       write(6,*)"xkapi_nco = ",xkapi_nco
c       write(6,*)istringer_nc,xnstari_nco,xnstare_nco
c       write(6,*)ztau_nco,zrho_nco,zfluxlim_nco
c      endif
c
c...   calculate neclassical conductivity (xkapi), diffusivity
c
       if(nkimod_nc.eq.4)then
         xkapi_nco=xkapi_nco*1.D-4*nim*1.D19  !1/(m*s)
         chiineo_m(jm)=xkapi_nco/(nim*1.D19)  !m**2/s
       else
c         xkapi_nco=xkapi_nco*1.D-4*nim*1.D19/elong_exp(jm)    !1/(m*s)
c         chiineo_m(jm)=xkapi_nco/nim/1.D19/gradrhosq_exp(jm)  !m**2/s
         xkapi_nco=xkapi_nco*1.D-4*nim*1.D19  !1/(m*s)
         chiineo_m(jm)=xkapi_nco/(nim*1.D19)  !m**2/s
       endif
       chieneo_m(jm)=chiineo_m(jm)/42.D0/DSQRT(amassgas_exp)
       deneo_m(jm)=chieneo_m(jm)
c compute toroidal momentum diffusivity
       lnlamda=15.94D0-0.5*LOG(nem)+LOG(tem)
c braginski ion collision time
       taui = 0.0661*(tim**1.5)/(nim*lnlamda)
       rhoi = rhosdam*a_unit_exp*SQRT(2.0*tim/tem)
       etaphineo_m(jm) = (3.0/5.0)*rhoi**2/taui
      elseif (use_xneo_m.eq.2) then
c       write(*,*) 'calling forcebal ...'
c       call forcebal
        deneo_m(jm)=deneo_exp(jm)
        dineo_m(jm)=dineo_exp(jm)
        chieneo_m(jm)=chieneo_exp(jm)
        chiineo_m(jm)=chiineo_exp(jm)
      endif
c
      chiineogb_m(jm)=chiineo_m(jm)/cgyrobohm_m(jm)
c
c  model for neoclassical diffusivities
c
      diffnem = xparam_pt(3)*deneo_m(jm)
      chietem = xparam_pt(4)*chieneo_m(jm)
      chiitim = xparam_pt(5)*chiineo_m(jm)
      etaphim = xparam_pt(6)*etaphineo_m(jm)
      if(imodel.eq.82)then
       diffnem = diffnem*drhodr(jm)
       chietem = chietem*drhodr(jm)
       chiitim = chiitim*drhodr(jm)
       etaphim = etaphim*drhodr(jm)
      endif
c
c     compute neoclassical flows
c
      if(ncl_flag.eq.0)then
c Large aspect ratio formulas from Y.B. Kim Phys. FLuids B3 (1991) 2050.
        cnc = -1.0*alpha_dia/bt_exp
        thetam=theta_exp(jm)
        rmajm = rmaj_exp(jm)
        rhom=arho_exp*rho(jm)
        rminm=rmin_exp(jm)
        qm=q_exp(jm)
        r_eps = rminm/rmajm
        fcm = 1.0+(0.46*r_eps-1.46)*DSQRT(r_eps)
        g_ft = (1.0 - fcm)/fcm
        kpol1=0.8839*fcm/(0.3477+0.4058*fcm)
        kpol2=(1.0 -fcm)/(1.0 + 1.167*fcm)
        kpol_m(jm)=0.8839*fcm/(0.3477+0.4058*fcm)
c calculate coefficient of poloidal collisional damping rate
        lnlamda=24.D0-DLOG(DSQRT(1.D13*nem)/(1000.0*tem))
        lnlamda=DMAX1(lnlamda,1.D0)
        taui=2.09D7*DSQRT(amassgas_exp)*
     >  ((tim*1000.0)**1.5)/(1.0D13*nim*lnlamda)  !s
        vthi=9.79D3*DSQRT(2000.0*tim/amassgas_exp)  !m/s
        wthi= vthi/(rmajm*qm)
        nu_star = 1.0/(wthi*taui*r_eps**1.5)
        zeffm = zeff_exp(jm)
        d_zeff = 2.23 + 5.32*zeffm+2.04*zeffm**2
        nu_norm = (nim*amassgas_exp*1.6726D-8)/taui  ! kg/m^3/s
c   normalized mu_00 from Table 1
c   banana regime    
        nu_banana=nu_norm*(0.53 + zeffm)
c   plateau regime  
        nu_plateau = nu_norm*3.54
c   PS regime 
        nu_braginski = nu_norm*(3.02 + 4.25*zeffm)/d_zeff
        if(ipert_gf.eq.0)then
c formula C21
          nu_pol_m(jm)= g_ft*nu_banana/
     > ((1+2.92*nu_star*nu_banana/nu_plateau)*
     > (1+nu_plateau/(6*wthi*taui*nu_braginski)))
c          write(*,*)jm,"nu_pol_m=",nu_pol_m(jm)
        endif
c calculate vneo and vdia on the half grid
c neoclassical poloidal velocity (km/sec)
        vneo(1) = -alpha_dia*curden_exp(jm)/(1.6022*nem*1000.0)
        vneo(2) = -(cnc/thetam)*kpol1*gradtim
        vneo(3) = -(cnc/thetam)*((kpol1+1.5*kpol2
     >     -1.0+1.0/zimp_exp)*gradtim - gradnim+ gradnzm/zimp_exp)
      endif
c
c diamagnetic velocity (km/sec)
c
      vdia(1) = (-cnc/thetam)*(gradtem+tem*gradnem/nem)
      vdia(2) = (cnc/thetam)
     >   *(gradtim +tim*gradfim/fim+tim*gradnem/nem)/zgas_exp
      vdia(3) = (cnc/thetam)*(gradtim+tim*gradfzm/fzm
     >   +tim*gradnem/nem)/zimp_exp
c
      if(jm.eq.1)then
        nu_pol_m(0)=nu_pol_m(1)
        kpol_m(0)=kpol_m(1)
      endif
c
c compute neoclassical fluxes
c
       neflux_neo = (-1.6022D-3)*diffnem*gradnem
       teflux_neo = (-1.6022D-3)*nem*chietem*gradtem
       tiflux_neo = (-1.6022D-3)*nim*chiitim*gradtim
c use large ExB rotation form
      ctorm = (c_tor(jm+1)+c_tor(jm))/2.0
      grad_c_tor=(c_tor(jm+1)-c_tor(jm))/dr(jm,2)
      vphiflux_neo=(-1.6726D-8)*amassgas_exp*nim*etaphim
     >     *cv*(ctorm*gradvexbm+grad_c_tor*vexbm)
      vparflux_neo = vphiflux_neo*c_per(jm)/c_tor(jm)
      niflux_neo = neflux_neo
      nzflux_neo = 0.0
      tzflux_neo = 0.0
      vparzflux_neo = 0.0
      vphizflux_neo = 0.0
      neflux_neo = neflux_neo*drhodr(jm)
      teflux_neo = teflux_neo*drhodr(jm)
      tiflux_neo = tiflux_neo*drhodr(jm)
      vphiflux_neo = vphiflux_neo*drhodr(jm)
c
c
      RETURN
      END !SUBROUTINE get_kapisn
c
      SUBROUTINE get_neo
      
      USE neo_interface
      IMPLICIT none
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
!
      real :: n0,a0,T0,v0,m0,cnc,thetam,cvpol     
!
      call xptor_neo_map
      call neo_run     
!
! neo normalizations
!
!      m0 = amassgas_exp                                !proton mass
!      n0 = nem                                         !10^19/m^3
!      a0 = rmin_exp(mxgrid)                            !m
!      T0 = tim                                         !kev
      m0=2.0
      n0 =1.0
      T0=1.0
      a0=1.0
      v0 = 9.79D3*DSQRT(T0*1.D3)/DSQRT(m0)   !m/sec
!     
      cnc = -1.0*alpha_dia/bt_exp
      thetam=theta_exp(jm)
      nu_pol_m(jm)=0.0
      if(ipert_gf.eq.0)then
        nu_pol_m(jm) = 1.6726D-8*neo_nclassvis_out(2)*m0*n0*v0/a0  ! nclass viscosity
c        write(*,*)jm,"nu_pol_m=",nu_pol_m(jm),neo_nclassvis_out(2)
      endif
      kpol_m(jm) = 0.0
!
! neoclassical poloidal velocity (km/sec)
! note that neo uses a coordinate system with theta and phi reversed from XPTOR
! here vneo is the toroidal rotation part <R*U_phi>/R0 = vphi = c_per*(vpol+vneo)+c_tor*(vexb+vdia)
! so vneo = sign_bt_exp*bt_exp*K and 
! neo_vpol = -ABS(Bp0)*sign_It_exp*K where the - sign comes fromt the reversed theta direction
! because vneo_pol = voloidal velocity at the outboard midplane in the direction of the neo theta.
!
      cvpol=-sign_It_exp*sign_Bt_exp*alpha_dia*(v0/cv)
     >       *ABS(bt_exp/Bp0(jm))
      vneo(1) = cvpol*neo_vpol_dke_out(1)
      vneo(2) = cvpol*neo_vpol_dke_out(2)
      if(neo_n_species_in.gt.2)vneo(3) = cvpol*neo_vpol_dke_out(3)
c
c diamagnetic velocity (km/sec)
c
      vdia(1) = (-cnc/thetam)*(gradtem+tem*gradnem/nem)
      vdia(2) = (cnc/thetam)
     >   *(gradtim +tim*gradfim/fim+tim*gradnem/nem)/zgas_exp
      vdia(3) = (cnc/thetam)*(gradtim+tim*gradfzm/fzm
     >   +tim*gradnem/nem)/zimp_exp
c
      if(jm.eq.1)then
        nu_pol_m(0)=nu_pol_m(1)
        kpol_m(0)=kpol_m(1)
      endif
c
      if(ineo.eq.-4.and.jm.le.jin_m)then
c skip neoclassical fluxes
        neflux_neo = 0.0
        niflux_neo = 0.0
        nzflux_neo = 0.0
        teflux_neo = 0.0
        tiflux_neo = 0.0
        tzflux_neo = 0.0
        vphiflux_neo = 0.0
        vphizflux_neo = 0.0
      else
c
c compute neoclassical fluxes
c
         neflux_neo = drhodr(jm)*(1.6022D-3)*n0*v0*
     >     (neo_pflux_dke_out(1)+neo_pflux_gv_out(1))
         niflux_neo = drhodr(jm)*(1.6022D-3)*n0*v0*
     >     (neo_pflux_dke_out(2)+neo_pflux_gv_out(2))
         teflux_neo = drhodr(jm)*(1.6022D-3)*n0*v0*T0*
     >     (neo_efluxtot_dke_out(1)+neo_efluxtot_gv_out(1))
         tiflux_neo = drhodr(jm)*(1.6022D-3)*n0*v0*T0*
     >     (neo_efluxtot_dke_out(2)+neo_efluxtot_gv_out(2))
         vphiflux_neo = -sign_It_exp*drhodr(jm)*(1.6726D-8)*
     > m0*n0*a0*v0*v0*(neo_mflux_dke_out(2)+neo_mflux_gv_out(2))
         if(neo_n_species_in.gt.2)then
           nzflux_neo = drhodr(jm)*(1.6022D-3)*n0*v0*
     >     (neo_pflux_dke_out(3)+neo_pflux_gv_out(3))
           tzflux_neo = drhodr(jm)*(1.6022D-3)*n0*v0*T0*
     >     (neo_efluxtot_dke_out(3)+neo_efluxtot_gv_out(3))
           vphizflux_neo = -sign_It_exp*drhodr(jm)*(1.6726D-8)*
     >   m0*n0*a0*v0*v0*(neo_mflux_dke_out(3)+neo_mflux_gv_out(3))
        endif
      endif
c      
      RETURN
      END ! SUBROUTINE get_neo
      SUBROUTINE  tglf_dv
c
      USE tglf_tg
      USE tglf_pkg
      IMPLICIT NONE
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      integer i,j,k,ik
      integer nk_e,nx_e
      real*8 exchgb_local,gfac,diffgb_local
      real*8 fc,vnewk3x,alpha_neo_hold,akappa1
      real*8 aiwt_jp1,xnimp_jp1,xnimp
      real*8 rhosdam,rmajm,rminm,rhom
      real*8 fac_imp_flow
      real*8 zeffm,qm, gb_unit,b_unit
      REAL*8 dtest(5,5),vtest(5,5)
      REAl*8 debyelorhos
      real*8 dr_loc,xnuei,taue,lnlamda
      real*8 vneom(nspecies),vdiam(nspecies)
      real*8 gradvneom(nspecies),gradvdiam(nspecies)
      real*8 gradvphim
      real*8 cparm,cperm,ctorm,grad_c_par,grad_c_per,grad_c_tor
      real*8 apolm,atorm,grad_a_pol,grad_a_tor
      real*8 vphim,mass_factor,vpartotm
      real*8 chi_DR,cnc,wstar0,wstar1
      real*8 dptotj1,dptotj0
c
c NOTE:
c The a_unit_exp convention is rmin at separatrix for TGLF for all geometries.
c This differes from GLF23 which uses a_unit_exp = arho_exp. This does not impact the 
c predicted plasma profiles since for both models the local flux does not depend
c on a_unit_exp so it is arbitrary
c
c      a_unit_exp = rmin_exp(mxgrid)
c      if(igeo_tg.eq.0)a_unit_exp=arho_exp
c
      cxnu = 1.0
      mass_factor = (1.0 + 
     >   amassimp_exp*fz_m(jm)/(amassgas_exp*fi_m(jm)))
      qm =(q_exp(jm+1)+q_exp(jm))/2.D0
      rmajm=(rmaj_exp(jm+1)+rmaj_exp(jm))/2.D0
      rhom=arho_exp*(rho(jm+1)+rho(jm))/2.D0
c recompute drhodr just to make sure
      drhodr(jm) = arho_exp*(rho(jm+1)-rho(jm))/
     >             (rmin_exp(jm+1)-rmin_exp(jm))
      rminm=(rmin_exp(jm+1)+rmin_exp(jm))/2.D0
      zeffm=(zeff_exp(jm+1)+zeff_exp(jm))/2.D0
c for TGLF derivatives are with respect to rmin for all geometries
      zpmne = -a_unit_exp*gradnem/nem
      zpmni = -a_unit_exp*gradnim/nim
      zpmnz = -a_unit_exp*gradnzm/nzm
      zpmte = -a_unit_exp*gradtem/tem
      zpmti = -a_unit_exp*gradtim/tim
      zpmne = drhodr(jm)*zpmne
      zpmni = drhodr(jm)*zpmni
      zpmnz = drhodr(jm)*zpmnz
      zpmte = drhodr(jm)*zpmte
      zpmti = drhodr(jm)*zpmti
c
      csdam=9.79D5*DSQRT(tem*1.D3)/
     > (a_unit_exp*100.D0)/DSQRT(amassgas_exp)
      csda_m(jm) = csdam
c      vdiam(1) = (vdia(1)+vdia_m(1,jm))/2.0
c      vdiam(2) = (vdia(2)+vdia_m(2,jm))/2.0
c      vdiam(3) = (vdia(3)+vdia_m(3,jm))/2.0
c      vneom(1) = (vneo(1)+vneo_m(1,jm))/2.0
c      vneom(2) = (vneo(2)+vneo_m(2,jm))/2.0
c      vneom(3) = (vneo(3)+vneo_m(3,jm))/2.0
c      gradvdiam(1) =(vdia(1)-vdia_m(1,jm))/dr(jm,2)
c      gradvdiam(2) =(vdia(2)-vdia_m(2,jm))/dr(jm,2)
c      gradvdiam(3) =(vdia(3)-vdia_m(3,jm))/dr(jm,2)
c      gradvneom(1) =(vneo(1)-vneo_m(1,jm))/dr(jm,2)
c      gradvneom(2) =(vneo(2)-vneo_m(2,jm))/dr(jm,2)
c      gradvneom(3) =(vneo(3)-vneo_m(3,jm))/dr(jm,2)
      vdiam(1) = (vdia_m(1,jm+1)+vdia_m(1,jm))/2.0
      vdiam(2) = (vdia_m(2,jm+1)+vdia_m(2,jm))/2.0
      vdiam(3) = (vdia_m(3,jm+1)+vdia_m(3,jm))/2.0
      vneom(1) = (vneo_m(1,jm+1)+vneo_m(1,jm))/2.0
      vneom(2) = (vneo_m(2,jm+1)+vneo_m(2,jm))/2.0
      vneom(3) = (vneo_m(3,jm+1)+vneo_m(3,jm))/2.0
      gradvdiam(1) =(vdia_m(1,jm+1)-vdia_m(1,jm))/dr(jm,2)
      gradvdiam(2) =(vdia_m(2,jm+1)-vdia_m(2,jm))/dr(jm,2)
      gradvdiam(3) =(vdia_m(3,jm+1)-vdia_m(3,jm))/dr(jm,2)
      gradvneom(1) =(vneo_m(1,jm+1)-vneo_m(1,jm))/dr(jm,2)
      gradvneom(2) =(vneo_m(2,jm+1)-vneo_m(2,jm))/dr(jm,2)
      gradvneom(3) =(vneo_m(3,jm+1)-vneo_m(3,jm))/dr(jm,2)
      cparm = (c_par(jm+1)+c_par(jm))/2.0
      cperm = (c_per(jm+1)+c_per(jm))/2.0
      ctorm = (c_tor(jm+1)+c_tor(jm))/2.0
      apolm = (a_pol(jm+1)+a_pol(jm))/2.0
      atorm = (a_tor(jm+1)+a_tor(jm))/2.0
      grad_c_par=(c_par(jm+1)-c_par(jm))/dr(jm,2)
      grad_c_per=(c_per(jm+1)-c_per(jm))/dr(jm,2)
      grad_c_tor=(c_tor(jm+1)-c_tor(jm))/dr(jm,2)
      grad_a_pol=(a_pol(jm+1)-a_pol(jm))/dr(jm,2)
      grad_a_tor=(a_tor(jm+1)-a_tor(jm))/dr(jm,2)
c
c      gamma_p_m(jm) = -(cv/(csdam))*drhodr(jm)*
c     >  (apolm*(gradvpolm+gradvneom(2))
c     >  +atorm*(gradvexbm+gradvdiam(2))
c     >  +grad_a_pol*(vpolm+vneom(2))
c     >  +grad_a_tor*(vexbm+vdiam(2)))
c      egamma_m(jm) = -cv/csdam*(rminm/rhom)*drhodr(jm)
c     >  *theta_exp(jm)*gradvexbm
c      egamma_m(jm) = -cv/csdam*(rminm/rhom)*drhodr(jm)
c     >  *theta_exp(jm)*(vexb_m(jm+1)-vexb_m(jm))/dr(jm,2)
c      if(jm.eq.ngrid-1.and.itport_pt(5).eq.0)then
c        egamma_m(jm)=0.0
c      endif
c       gradvphim = -gamma_p_m(jm)*csdam/cv
c       gradvphim = (cv/csdam)*drhodr(jm)*gradvexbm
c       gradvphim =
c     >  cperm*(gradvpolm+gradvneom(2))+ctorm*(gradvexbm+gradvdiam(2))
c     > +grad_c_per*(vpolm+vneom(2))+grad_c_tor*(vexbm+vdiam(2))
c
c   local magnetic field unit 
      if(igeo_tg.eq.0)then
        b_unit = bt_exp
c        if(bt_flag.gt.0)b_unit = bteff_exp(jm)
      else
        b_unit = bt_exp*(rhom/rminm)*drhodr(jm)
      endif
c   local rho_star
      rhosdam=(1.02D2*DSQRT((tem*1.D3)*amassgas_exp)/
     >  (b_unit*1.D4))/(a_unit_exp*100.D0)
      rhosda_m(jm) = rhosdam
c   local gyrobohm unit of diffusion
      cgyrobohm_m(jm)= csdam*(rhosdam*a_unit_exp)**2
c   local beta
      betae_m(jm) = 4.03D-3*nem*tem/(b_unit**2)
      if(iparam_pt(12).eq.1)then
        betae_m(jm) = 4.03D-3*ne_exp(jm)*te_exp(jm)/(b_unit**2)
      endif
c   debye length/rhos
      debyelorhos = 7.43D2*SQRT((tem*1.D3)/(nem*1.D13))/
     > (a_unit_exp*100.0*rhosdam)

c   gks collisionality nu_ei*a/vth_i  where vth_i = sqrt(2*ti/mi)
       vnewk3x= 0.5D0*(a_unit_exp)*DSQRT(amassgas_exp/2.D0)*
     >   0.117D0*nem/DSQRT(tim*tem**3)
c   switch to Cs/a units
       xnu_m(jm) =vnewk3x*DSQRT(2.D0*tim/tem)
c
c  new calculation of xnu_m with lnlamda 
c  xnuei = 3/4 sqt(Pi)/taue
c  lnlamda and taue from NRL formulary
c  note for Te=1Kev and ne=10**13 lnlamda = 15.94 and taue= 1.088D-3/lnlamda
       lnlamda=15.94D0-0.5*LOG(nem)+LOG(tem)
       taue = 1.088D-3*(tem**1.5)/(nem*lnlamda)
c  xnuei = 3/4 (Pi**0.5)/taue
       xnuei = 1.329/taue
c       if(jm.eq.35)write(*,*)"lnlamda=",lnlamda,
c     > "xnugks=",xnu_m(jm),"xnutglf=",xnuei/csdam
       xnu_m(jm) = xnuei/csdam
c
       if(ipert_gf.eq.0)alpha_m(jm)=drhodr(jm)*qm**2*rmajm
     >   *betae_m(jm)*((tim*nim/tem/nem)*
     >   (zpmni+zpmti)+zpmne+zpmte)/a_unit_exp
c bananna regime ie collisionless limit formulas
        alpha_neo_hold=alpha_neo
        fc=1-1.46D0*DSQRT(rminm/rmajm)+
     >      0.46D0*DSQRT(rminm/rmajm)**3
        akappa1=0.8839D0*fc/(0.3477D0+0.4058D0*fc)
        if(irot1.eq.1) alpha_neo=-akappa1+1.D0
c set TGLF gradients
      rlts_tg(1)=zpmte
      rlts_tg(2)=zpmti
      rlts_tg(3)=zpmti
      rlns_tg(1)=zpmne
      rlns_tg(2)=zpmni
      rlns_tg(3)=zpmnz
      vexb_shear_tg = -(cv/csdam)*(rminm/rhom)*drhodr(jm)
     >  *theta_exp(jm)*gradvexbm
      if(ipert_gf.eq.0)egamma_m(jm)=vexb_shear_tg*csdam/cv
      if(iexb.ge.1)then
c        vexb_shear_tg = egamma_exp(jm)
        vexb_shear_tg =(cv/csdam)* egamma_exp(jm)
      endif
      if(doppler_shear_model.eq.1)then
c this takes out the pertubation of the Doppler shear through gradvexbm
c        vexb_shear_tg = -(cv/csdam)*(rminm/rhom)*drhodr(jm)
c     >  *theta_exp(jm)*(vexb_m(jm+1)-vexb_m(jm))/dr(jm,2)
        vexb_shear_tg = (cv/csdam)*doppler_shear_m(jm)
      endif
c      if(jm.eq.ngrid-1.and.itport_pt(5).eq.0)vexb_shear_tg=0.0
c
c      vpar_shear_tg(1) = -sign_Bt_exp*(cv/(csdam))*drhodr(jm)*
c     >  (apolm*(gradvpolm+gradvneom(1))+atorm*(gradvexbm+gradvdiam(1))) 
cc     >  +(vpolm+vneom(1))*grad_a_pol+(vexbm+vdiam(1))*grad_a_tor)
c      vpar_shear_tg(2) = -sign_Bt_exp*(cv/(csdam))*drhodr(jm)*
c     >  (apolm*(gradvpolm+gradvneom(2))+atorm*(gradvexbm+gradvdiam(2))) 
cc     >  +(vpolm+vneom(2))*grad_a_pol+(vexbm+vdiam(2))*grad_a_tor)
c      vpar_shear_tg(3) = -sign_Bt_exp*(cv/(csdam))*drhodr(jm)*
c     >  (apolm*(gradvpolm+gradvneom(3))+atorm*(gradvexbm+gradvdiam(3))) 
cc     >  +(vpolm+vneom(3))*grad_a_pol+(vexbm+vdiam(3))*grad_a_tor)
      vpar_shear_tg(1) = -(cv/csdam)*drhodr(jm)*
     >  (rmajm/rmajor_exp)*(gradvexbm+gradvdiam(1)+
     >  (apolm/atorm)*(gradvpolm+gradvneom(1))) 
      vpar_shear_tg(2) = -(cv/csdam)*drhodr(jm)*
     >  (rmajm/rmajor_exp)*(gradvexbm+gradvdiam(2)+
     >  (apolm/atorm)*(gradvpolm+gradvneom(2))) 
      vpar_shear_tg(3) = -(cv/csdam)*drhodr(jm)*
     >  (rmajm/rmajor_exp)*(gradvexbm+gradvdiam(3)+
     >  (apolm/atorm)*(gradvpolm+gradvneom(3))) 
c
      if(ipert_gf.eq.0)gamma_p_m(jm)=vpar_shear_tg(2)
c
      cnc = -alpha_dia/bt_exp
      wstar0 = -(cnc/theta_exp(jm))*(te_m(jm+1)-te_m(jm-1))
     >   /((dr(jm,1)+dr(jm,2))*(te_m(jm+1)+te_m(jm-1))/2.0)
      wstar1 = -(cnc/theta_exp(jm+1))*(te_m(jm+2)-te_m(jm))
     >  /((dr(jm+1,1)+dr(jm+1,2))*(te_m(jm+2)+te_m(jm))/2.0)
      vts_shear_tg(1) =  -cv/csdam*(rminm/rhom)*drhodr(jm)
     >  *theta_exp(jm)*te_m(jm)*(wstar1 - wstar0)/dr(jm,2)
      wstar0 = -(cnc/theta_exp(jm))*(ti_m(jm+1)-ti_m(jm-1))
     >  /((dr(jm,1)+dr(jm,2))*(ti_m(jm+1)+ti_m(jm-1))/2.0)
      wstar1 = -(cnc/theta_exp(jm+1))*(ti_m(jm+2)-ti_m(jm))
     > /((dr(jm+1,1)+dr(jm+1,2))*(ti_m(jm+2)+ti_m(jm))/2.0)
      vts_shear_tg(2) =  -cv/csdam*(rminm/rhom)*drhodr(jm)
     >  *theta_exp(jm)*te_m(jm)*(wstar1 - wstar0)/dr(jm,2)
      vts_shear_tg(3) = vts_shear_tg(2)
      wstar0 = -(cnc/theta_exp(jm))*(ne_m(jm+1)-ne_m(jm-1))
     >  /((dr(jm,1)+dr(jm,2))*(ne_m(jm+1)+ne_m(jm-1))/2.0)
      wstar1 = -(cnc/theta_exp(jm+1))*(ne_m(jm+2)-ne_m(jm))
     > /((dr(jm+1,1)+dr(jm+1,2))*(ne_m(jm+2)+ne_m(jm))/2.0)
      vns_shear_tg(1) =  -cv/csdam*(rminm/rhom)*drhodr(jm)
     >  *theta_exp(jm)*te_m(jm)*(wstar1 - wstar0)/dr(jm,2) 
      wstar0 = -(cnc/theta_exp(jm))*(ni_m(jm+1)-ni_m(jm-1))
     >  /((dr(jm,1)+dr(jm,2))*(ni_m(jm+1)+ni_m(jm-1))/2.0)
      wstar1 = -(cnc/theta_exp(jm+1))*(ni_m(jm+2)-ni_m(jm))
     > /((dr(jm+1,1)+dr(jm+1,2))*(ni_m(jm+2)+ni_m(jm))/2.0)
      vns_shear_tg(2) =  -cv/csdam*(rminm/rhom)*drhodr(jm)
     >  *theta_exp(jm)*te_m(jm)*(wstar1 - wstar0)/dr(jm,2)
      wstar0 = -(cnc/theta_exp(jm))*(nz_m(jm+1)-nz_m(jm-1))
     >  /((dr(jm,1)+dr(jm,2))*(nz_m(jm+1)+nz_m(jm-1))/2.0)
      wstar1 = -(cnc/theta_exp(jm+1))*(nz_m(jm+2)-nz_m(jm))
     > /((dr(jm+1,1)+dr(jm+1,2))*(nz_m(jm+2)+nz_m(jm))/2.0)
      vns_shear_tg(3) =  -cv/csdam*(rminm/rhom)*drhodr(jm)
     >  *theta_exp(jm)*te_m(jm)*(wstar1 - wstar0)/dr(jm,2)
c
c      if(jm.eq.ngrid-1)then
c set diamagnetic and neoclassical flow gradeints to zero at boundary
c        vpar_shear_tg(1) = -sign_Bt_exp*(cv/(csdam))*drhodr(jm)*
c     >  (apolm*(gradvpolm)+atorm*(gradvexbm)) 
c     >  +(vpolm)*grad_a_pol+(vexbm)*grad_a_tor)
c        vpar_shear_tg(2)=vpar_shear_tg(1)
c        vpar_shear_tg(3)=vpar_shear_tg(1)
c        vns_shear_tg(1) = 0.0
c        vns_shear_tg(2) = 0.0
c        vns_shear_tg(3) = 0.0
c        vts_shear_tg(1) = 0.0
c        vts_shear_tg(2) = 0.0
c        vts_shear_tg(3) = 0.0
c      endif
      if(alpha_p_tg.eq.0.0)then
        vpar_shear_tg(1)=0.0
        vpar_shear_tg(2)=0.0
        vpar_shear_tg(3)=0.0
      endif
c initialized fluxes
      neflux_glf = 0.0
      niflux_glf = 0.0
      nzflux_glf = 0.0
      teflux_glf = 0.0
      tiflux_glf = 0.0
      tzflux_glf = 0.0
      vphiflux_glf = 0.0
      vphizflux_glf = 0.0
      vparflux_glf = 0.0
      vparzflux_glf = 0.0
c
      if(ineo.eq.-2.and.q_exp(jm).lt.1.0)go to 86 !skip TGLF
      if(ineo.eq.-3.and.interchange_DR_m(jm).lt.0.0)go to 86 ! skip TGLF
      if(ineo.eq.-4.and.jm.le.jin_m)go to 86
c
c  do not call TGLF to compute variations if the first call has no flux
c       if(ipert_gf.eq.0)write(*,*)"glf2d_dv",jm
cc       if(ipert_gf.ne.0)then
cc        if(v2_bar(jm).eq.0.0)go to 86  !skip TGLF fluxes
cc       endif

c  transfer inputs to TGLF
c
c
      CALL put_gaussian_width(width_max_tg,width_min_tg,nwidth_tg, 
     > find_width_tg)
c
c      CALL put_switches(iflux_tg,use_bper_tg,use_bpar_tg,
c     > use_mhd_rule_tg,use_bisection_tg, 
c     > ibranch_tg,nmodes_tg,nb_max_tg,nb_min_tg,nx_tg,nky_tg)
c
      CALL put_gradients(rlns_tg,rlts_tg,vpar_shear_tg,vexb_shear_tg)
c
c      CALL put_profile_shear(vns_shear_tg,vts_shear_tg)
c      if(ipert_gf.eq.0)then
c      write(*,*)"pro",jm,vns_shear_tg(2),vts_shear_tg(2),vexb_shear_tg
c      endif
c
      taus_tg(1)=1.0
      taus_tg(2)=tim/tem
      taus_tg(3)=tim/tem
      as_tg(1)=1.0
      as_tg(2)= nim/nem
      as_tg(3)= nzm/nem
      betae_tg=cbetae*betae_m(jm)
      xnue_tg =xnu_m(jm)
      zeff_tg=zeffm
      vexb_tg=0.0
c      vpar_tg(1) = sign_Bt_exp*cv*(a_pol(jm)*(vpolm+vneom(1))
c     >     +a_tor(jm)*(vexbm+vdiam(1)))/(a_unit_exp*csdam)
c      vpar_tg(2) = sign_Bt_exp*cv*(a_pol(jm)*(vpolm+vneom(2))
c     >     +a_tor(jm)*(vexbm+vdiam(2)))/(a_unit_exp*csdam)
c      vpar_tg(3) = sign_Bt_exp*cv*(a_pol(jm)*(vpolm+vneom(3))
c     >     +a_tor(jm)*(vexbm+vdiam(3)))/(a_unit_exp*csdam)
      vpar_tg(1) = cv*rmajm/rmajor_exp*
     >   (vexbm+vdiam(1)+(apolm/atorm)*(vpolm+vneom(1)))
     >   /(a_unit_exp*csdam)
      vpar_tg(2) = cv*rmajm/rmajor_exp*
     >   (vexbm+vdiam(2)+(apolm/atorm)*(vpolm+vneom(2)))
     >   /(a_unit_exp*csdam)
      vpar_tg(3) = cv*rmajm/rmajor_exp*
     >  (vexbm+vdiam(3)+(apolm/atorm)*(vpolm+vneom(3)))
     >   /(a_unit_exp*csdam)
      if(alpha_p_tg.eq.0.0)then
        vpar_tg(1)=0.0
        vpar_tg(2)=0.0
        vpar_tg(3)=0.0
        vexb_tg=0.0
      endif
c      debye_tg = debyelorhos**2
       debye_tg = debyelorhos
c      if(ipert_gf.eq.0)write(*,*)"debye_tg=",debye_tg
c      write(*,*)"jm=",jm,"ipert_gf=",ipert_gf
c      write(*,*)"debug","gradtem=",gradtem,"tem=",tem
c      write(*,*)"debug","gradtim=",gradtim,"tim=",tim     
c
      CALL put_averages(taus_tg,as_tg,vpar_tg,vexb_tg,betae_tg,xnue_tg,
     > zeff_tg,debye_tg)
c
c      if(ipert_gf.eq.0)then
       if(igeo_tg.eq.0)then
c s-alpha geometry
        rmin_tg = rminm/a_unit_exp
        rmaj_tg = rmajm/a_unit_exp
        q_tg = qm
        shat_tg = shat_exp(jm)
        alpha_tg=xalpha*alpha_exp(jm)
        if(ialphastab.gt.0) then
          alpha_tg=xalpha*alpha_m(jm)
        endif
        CALL put_s_alpha_geometry(rmin_tg,rmaj_tg,q_tg,shat_tg,alpha_tg, 
     >    xwell_tg,theta0_tg,b_model_tg,ft_model_tg)
c
       elseif(igeo_tg.eq.1)then  ! miller geometry
c
         dr_loc = (rmin_exp(jm+1)-rmin_exp(jm))/a_unit_exp
         rmin_tg = rminm/a_unit_exp
         rmaj_tg = rmajm/a_unit_exp
         zmaj_tg = 0.0
         q_tg = qm
         q_prime_tg = (q_tg/rmin_tg)*(q_exp(jm+1)-q_exp(jm))/dr_loc
c
         p_prime_tg = (q_tg/rmin_tg)*(1.6022D-4/b_unit**2)*
     >                (ptot_exp(jm+1)-ptot_exp(jm))/dr_loc
         if(ialphastab.eq.1)then
           p_prime_tg = (q_tg/rmin_tg)*(1.6022D-4/b_unit**2)*
     > (nem*gradtem+tem*gradnem+(nim+nzm)*gradtim+tim*(gradnim+gradnzm))
         dptotj1 = ptot_exp(jm+1)-(ne_exp(jm+1)*te_exp(jm+1)+
     >   (ni_exp(jm+1)+nz_exp(jm+1))*ti_exp(jm+1))
         dptotj0 = ptot_exp(jm)-(ne_exp(jm)*te_exp(jm)+
     >   (ni_exp(jm)+nz_exp(jm))*ti_exp(jm))
         p_prime_tg = p_prime_tg 
     >   +(q_tg/rmin_tg)*(1.6022D-4/b_unit**2)*(dptotj1-dptotj0)/dr_loc           
         endif
         p_prime_tg = xalpha*p_prime_tg
c
         drmindx_tg=1.0
         drmajdx_tg = (rmaj_exp(jm+1)-rmaj_exp(jm))/(dr_loc*a_unit_exp)
         dzmajdx_tg=0.0
         kappa_tg = 0.5*(elong_exp(jm+1)+elong_exp(jm))
         s_kappa_tg = (rmin_tg/kappa_tg)*
     >                (elong_exp(jm+1)-elong_exp(jm))/dr_loc
c         s_kappa_tg = dmax1(s_kappa_tg,1.e-2)
c why prevent negative s_kappa? 
         delta_tg = 0.5*(delta_exp(jm+1)+delta_exp(jm))
         s_delta_tg = rmin_tg*(delta_exp(jm+1)-delta_exp(jm))/dr_loc
         zeta_tg=0.0
         s_zeta_tg=0.0
        CALL put_Miller_geometry(rmin_tg,rmaj_tg,zmaj_tg,drmindx_tg,
     >  drmajdx_tg,dzmajdx_tg,kappa_tg,s_kappa_tg,delta_tg,s_delta_tg,
     >  zeta_tg,s_zeta_tg,q_tg,q_prime_tg,p_prime_tg,kx0_tg)
c
       else
        write(*,*)"igeo_tg invalid",igeo_tg
        stop
       endif
c      endif
c
c
c      write(*,*)"glf2d_dv",jm,"ipert_gf=",ipert_gf
c      if(ipert_gf.ne.0)new_eikonal_tg=.FALSE.
c use eikonal from last call to tglf_TM 
c      CALL put_eikonal(new_eikonal_tg)
      CALL tglf_TM
c      new_eikonal_tg=.TRUE.
c
      if(ipert_gf.eq.0.and.jm.eq.j_write_state)then
        call write_tglf_input      
      endif
      if(ipert_gf.eq.0)then
        v2_bar(jm) = get_v_bar_sum()
c        write(*,*)jm,"v2_bar=",v2_bar(jm)
        anrate_m(jm)= get_growthrate(1)
        anfreq_m(jm)= get_frequency(1)
c      write(*,*)"debug gamma",anrate_m(jm)
c      write(*,*)"debug freq",anfreq_m(jm)
c      write(*,*)"debug R_unit=",get_R_unit(),"q_unit=",get_q_unit()
      endif
c
c gb_unit = common gyrobohm factor for fluxes = cs*(rhos/a)**2
      gb_unit = cgyrobohm_m(jm)/a_unit_exp
cgms      if(igeo_tg.ne.0)gb_unit = gb_unit*drhodr(jm)
      neflux_glf = 1.6022D-3*nem*gb_unit*get_particle_flux(1,1)
      niflux_glf = 1.6022D-3*nem*gb_unit*get_particle_flux(2,1)
      teflux_glf = 1.6022D-3*nem*tem*gb_unit*get_energy_flux(1,1)
      tiflux_glf = 1.6022D-3*nem*tem*gb_unit*get_energy_flux(2,1)
      vparflux_glf = (1.6726D-8)  
     >   *nem*gb_unit*csdam*a_unit_exp*get_stress_par(2,1)
     >   *amassgas_exp/(ABS(bt_exp/B_unit))
      vphiflux_glf = (1.6726D-8)*amassgas_exp*a_unit_exp
     >  *nem*gb_unit*csdam*a_unit_exp*get_stress_tor(2,1)
      exch_glf(jm)=1.6022D-3*nem*tem*gb_unit*get_exchange(1,1)
     >             /a_unit_exp
      if(ns_tg.eq.3)then
        nzflux_glf = 1.6022D-3*nem*gb_unit*get_particle_flux(3,1)
        tzflux_glf = 1.6022D-3*nem*tem*gb_unit*get_energy_flux(3,1)
        vphizflux_glf = (1.6726D-8)*amassgas_exp*a_unit_exp
     >  *nem*gb_unit*csdam*a_unit_exp*get_stress_tor(3,1)  
        vparzflux_glf = (1.6726D-8) 
     >   *nem*gb_unit*csdam*a_unit_exp*get_stress_par(3,1)
     >   *amassgas_exp/(ABS(bt_exp/B_unit))
      endif
      if(use_bper_tg)then
        neflux_glf = neflux_glf 
     >    + 1.6022D-3*nem*gb_unit*get_particle_flux(1,2)
        teflux_glf = teflux_glf 
     >    + 1.6022D-3*nem*tem*gb_unit*get_energy_flux(1,2)
        tiflux_glf = tiflux_glf 
     >    + 1.6022D-3*nem*tem*gb_unit*get_energy_flux(2,2)
        vparflux_glf = vparflux_glf + (1.6726D-8)  
     >   *nem*gb_unit*csdam*a_unit_exp*get_stress_par(2,2)
     >   *amassgas_exp/(ABS(bt_exp/B_unit))
        vphiflux_glf = vphiflux_glf 
     >   + (1.6726D-8)*amassgas_exp*a_unit_exp
     >    *nem*gb_unit*csdam*a_unit_exp*get_stress_tor(2,2)
        if(ns_tg.eq.3)then
          tiflux_glf = tiflux_glf 
     >    +1.6022D-3*nem*tem*gb_unit*get_energy_flux(3,2)
          vphizflux_glf = vphizflux_glf 
     >     +(1.6726D-8)*amassgas_exp*a_unit_exp
     >    *nem*gb_unit*csdam*a_unit_exp*get_stress_tor(3,2)  
          vparzflux_glf = vparzflux_glf +(1.6726D-8) 
     >    *nem*gb_unit*csdam*a_unit_exp*get_stress_par(3,2)
     >    *amassgas_exp/(ABS(bt_exp/B_unit))
        endif
      endif
      if(use_bpar_tg)then
        neflux_glf = neflux_glf 
     >   +  1.6022D-3*nem*gb_unit*get_particle_flux(1,3)
        teflux_glf = teflux_glf 
     >   +  1.6022D-3*nem*tem*gb_unit*get_energy_flux(1,3)
        tiflux_glf = tiflux_glf 
     >   +  1.6022D-3*nem*tem*gb_unit*get_energy_flux(2,3)
        vparflux_glf = vparflux_glf + (1.6726D-8)  
     >   *nem*gb_unit*csdam*a_unit_exp*get_stress_par(2,3)
     >   *amassgas_exp/(ABS(bt_exp/B_unit))
        vphiflux_glf = vphiflux_glf 
     >   + (1.6726D-8)*amassgas_exp*a_unit_exp
     >    *nem*gb_unit*csdam*a_unit_exp*get_stress_tor(2,3)
        if(ns_tg.eq.3)then
          tiflux_glf = tiflux_glf 
     >    + 1.6022D-3*nem*tem*gb_unit*get_energy_flux(3,3)
          vphizflux_glf = vphizflux_glf 
     >     +(1.6726D-8)*amassgas_exp*a_unit_exp
     >    *nem*gb_unit*csdam*a_unit_exp*get_stress_tor(3,3)  
          vparzflux_glf = vparzflux_glf +(1.6726D-8) 
     >    *nem*gb_unit*csdam*a_unit_exp*get_stress_par(3,3)
     >    *amassgas_exp/(ABS(bt_exp/B_unit))
        endif
      endif
       tefluxm_etg = 1.6022D-3*nem*tem*gb_unit*(get_energy_flux(1,1)
     >              +get_energy_flux(1,2) -get_q_low(1))
       if(ns_tg.eq.2)then
c model for impurity contributions to viscous stress
         vphiflux_glf = mass_factor*vphiflux_glf
         vparflux_glf = mass_factor*vparflux_glf
         vphizflux_glf = 0.0
         vparzflux_glf = 0.0
         niflux_glf = neflux_glf
         nzflux_glf = 0.0
c model for impurity contribution to ion energy flux
         tiflux_glf = (1.0+nzm/nim)*tiflux_glf
         tzflux_glf = 0.0
       endif
c
c      write(*,*)jm,"ipert = ",ipert_gf,"v2_bar =",v2_bar(jm)
c      write(*,*)"debug Qion",get_energy_flux(2)
c      write(*,*)"debug teflux_glf=",teflux_glf
c
 86   continue   !skip TGLF fluxes
c
c
c multiply fluxes by drhdr factor from volume derivative dV/dr=dV/drho*drhodr
c to map from r-grid fluxes to rho-grid fluxes
c
      neflux_glf = neflux_glf*drhodr(jm)
      niflux_glf = niflux_glf*drhodr(jm)
      nzflux_glf = nzflux_glf*drhodr(jm)
      teflux_glf = teflux_glf*drhodr(jm)
      tiflux_glf = tiflux_glf*drhodr(jm)
      tzflux_glf = tzflux_glf*drhodr(jm)
      vphiflux_glf = vphiflux_glf*drhodr(jm)
      vphizflux_glf = vphizflux_glf*drhodr(jm)
      vparflux_glf = vparflux_glf*drhodr(jm)
      vparzflux_glf = vparzflux_glf*drhodr(jm)
c
      if(i_proc.eq.0.and.jm.eq.-1)then
       open(unit=11,file="gyro_input",status="unknown")
       write(11,*)"gyro input for shot ",shot
       write(11,*)"RADIUS = ",rmin_tg
       write(11,*)"ASPECT_RATIO = ",rmaj_tg
       write(11,*)"KAPPA = ",kappa_tg
       write(11,*)"S_KAPPA = ",s_kappa_tg
       write(11,*)"DELTA0 = ",delta_tg
       write(11,*)"DELTA1 = ",0.0
       write(11,*)"S_DELTA0 = ",s_delta_tg*SQRT(1.0-delta_tg**2) 
       write(11,*)"S_DELTA1 = ",0.0
       write(11,*)"SHIFT = ",drmajdx_tg
       write(11,*)"SAFETY_FACTOR = ",q_tg
       write(11,*)"SHEAR = ",q_prime_tg*(rmin_tg/q_tg)**2
       write(11,*)"RHO_STAR = ",rhosda_m(jm)
       write(11,*)"Z_EFF = ",zeff_tg
       write(11,*)"MACH = ",vpar_tg(2)
       write(11,*)"vpar_tg(1)=",vpar_tg(1)
       write(11,*)"vpar_tg(3)=",vpar_tg(3)
       write(11,*)"PGAMMA = ",vpar_shear_tg(2)
       write(11,*)"vpar_shear_tg(1)=",vpar_shear_tg(1)
       write(11,*)"vpar_shear_tg(3)=",vpar_shear_tg(3)
       write(11,*)"GAMMA_E = ",vexb_shear_tg
       write(11,*)"LAMDA_DEBYE = ",debye_tg
       write(11,*)"NU_EI = ",xnue_tg
       write(11,*)"NU_E_DRAG = ",0.0
       write(11,*)"NU_I_KROOK = ",0.0
       write(11,*)"# ION 1"
       write(11,*)"NI_OVER_NE = ",as_tg(2)/as_tg(1)
       write(11,*)"TI_OVER_TE = ",taus_tg(2)/taus_tg(1)
       write(11,*)"DLNNDR = ",rlns_tg(2)
       write(11,*)"DLNTDR = ",rlts_tg(2)
       write(11,*)"Z = ",zs_tg(2)
       write(11,*)"Electrons"
       write(11,*)"DLNNDR_ELECTRON = ",rlns_tg(1)
       write(11,*)"DLNTDR_ELECTRON = ",rlts_tg(1)
       write(11,*)"Z_ELECTRON = ",zs_tg(1)
       write(11,*)"BETAE_UNIT = ",betae_tg
       write(11,*)"MU_ELECTRON = ",SQRT(mass_tg(2)/mass_tg(1))
       write(11,*)"---- EFFECTIVE PHYSICAL UNITS ----"
       write(11,*)"b_unit_norm = ",b_unit
       write(11,*)"a_unit_norm = ",a_unit_exp
       write(11,*)"csda_norm = ",csdam
       write(11,*)"cs_norm = ",csdam*a_unit_exp
       write(11,*)"tem_norm = ",tem
       write(11,*)"den_norm = ",nem
       write(11,*)"chi_gb_norm = ",cgyrobohm_m(jm)
       write(11,*)"bt_exp = ",bt_exp
       write(11,*)"ion #1 mass/mp = ",amassgas_gf
       write(11,*)"TGLF OUTPUT"
       write(11,*)"particle_flux = ",get_particle_flux(1,1)
       write(11,*)"electron energy flux = ",get_energy_flux(1,1)
       write(11,*)"ion energy flux = ",get_energy_flux(2,1)
       write(11,*)"toroidal stress = ",get_stress_tor(2,1)
       close(11)
      endif
c
c
c
c      write(*,*)"debug tglf",jm
c
c      write(*,*) 'particle_flux(1) = ',get_particle_flux(1)
c      write(*,*) 'particle_flux(2) = ',get_particle_flux(2)
c      write(*,*) 'energy_flux(1) = ',get_energy_flux(1)
c      write(*,*) 'energy_flux(2) = ',get_energy_flux(2)
c      write(*,*) 'n_bar_sum(1) = ',get_n_bar_sum(1)
c      write(*,*) 'n_bar_sum(2) = ',get_n_bar_sum(2)
c      write(*,*) 't_bar_sum(1) = ',get_t_bar_sum(1)
c      write(*,*) 't_bar_sum(2) = ',get_t_bar_sum(2)
c      write(*,*) 'phi_bar_sum = ',get_phi_bar_sum()
c
c  
      END   !SUBROUTINE tglf_dv
c
      SUBROUTINE ADHOC_DV
c
      IMPLICIT NONE
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      real*8 art_diff
      real*8 ctorm,grad_c_tor,gradvphim
      real*8 diff_adhoc,chie_adhoc,chii_adhoc,eta_tor_adhoc
c
      diffnem = 0.0
      chietem = 0.0
      chiitim = 0.0
      etaphim = 0.0
      diff_adhoc = ABS(diff_exp(jm))
      chie_adhoc = ABS(chie_exp(jm))
      chii_adhoc = ABS(chii_exp(jm))
      eta_tor_adhoc = ABS(eta_tor_exp(jm))
c
c add in artificial diffusion terms
c
      if(adiff_dv.gt.0.0)then
       art_diff=adiff_dv*cgyrobohm_m(jm)
       diffnem = diffnem 
     > + art_diff*(dr(jm,2)*DABS(zpmne)/arho_exp)**2
       chietem = chietem 
     > + art_diff*DABS(dr(jm,2)*zpmte/arho_exp)**2
       chiitim = chiitim 
     > + art_diff*DABS(dr(jm,2)*zpmti/arho_exp)**2
       etaphim = etaphim 
     > + 1.D-2*art_diff*DABS(gamma_p_m(jm)/30.0)**2
      endif
c
c sawtooth region enhanced neoclassical model
      if((ineo.eq.-1.or.ineo.le.-2).and.q_exp(jm).lt.1.0)then
          diffnem = diffnem + chiineo_m(jm)
          if(imodel.ne.81)chietem = chietem + chiineo_m(jm)
c          chiitim = chiitim + chiineo_m(jm)
          etaphim = etaphim + chiineo_m(jm)
      endif
c
      if(ineo.eq.-3.and.interchange_DR_m(jm).lt.0.0)then
c use simple model for interchange mode fluxes
        diffnem = diffnem + diff_adhoc*drhodr(jm)*drhodr(jm)
        chietem = chietem + chie_adhoc*drhodr(jm)*drhodr(jm)
        chiitim = chiitim + chii_adhoc*drhodr(jm)*drhodr(jm)
        etaphim = etaphim + eta_tor_adhoc*drhodr(jm)*drhodr(jm)
      endif
c
      if(ineo.eq.-4.and.jm.le.jin_m)then
c use experimental diffusivities near the magnetic axis
        diffnem = diffnem + diff_adhoc*drhodr(jm)*drhodr(jm)
        chietem = chietem + chie_adhoc*drhodr(jm)*drhodr(jm)
        chiitim = chiitim + chii_adhoc*drhodr(jm)*drhodr(jm)
        etaphim = etaphim + eta_tor_adhoc*drhodr(jm)*drhodr(jm)
      endif
c
        neflux_adhoc = -1.6022D-3*diffnem*gradnem
        teflux_adhoc = -1.6022D-3*nem*chietem*gradtem
        tiflux_adhoc = -1.6022D-3*nim*chiitim*gradtim
        ctorm = (c_tor(jm+1)+c_tor(jm))/2.0
        grad_c_tor=(c_tor(jm+1)-c_tor(jm))/dr(jm,2)
c        gradvphim =cv*(ctorm*gradvexbm+grad_c_tor*vexbm)
        gradvphim =cv*ctorm*gradvexbm
        vphiflux_adhoc = -1.6726D-8*amassgas_exp*nim*
     >                    rmajor_exp*etaphim*gradvphim
        vparflux_adhoc = vphiflux_adhoc*c_per(jm)/c_tor(jm)
        niflux_adhoc = neflux_adhoc
        nzflux_adhoc = 0.0
        tzflux_adhoc = 0.0
        vparzflux_adhoc = 0.0
        vphizflux_adhoc = 0.0
c
c
      RETURN
      END  !ADHOC_DV
c
      SUBROUTINE test_dv
c
c test for debug of d,v
c
      IMPLICIT NONE
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/share.m'
      include '../inc/sharegk.m'
      include '../inc/ptor.m'
      include '../inc/glf.m'
c
      integer i,j,k
      REAL*8 dtest(5,5),vtest(5,5)
c
      do i=1,5
      do j=1,5
         dtest(i,j) = REAL(i+j)
         vtest(i,j) = 2.0*dtest(i,j)
      enddo
      enddo
      nefluxm=0.0
      tefluxm=0.0
      tifluxm=0.0
      vphifluxm=0.0
      vparfluxm=0.0
      i=0
      if(itport_pt(1).ne.0)then
        i=i+1
        j=0
        if(itport_pt(1).ne.0)then
          j=j+1
          nefluxm = nefluxm -dtest(i,j)*gradnem + vtest(i,j)*nem
        endif
        if(itport_pt(2).ne.0)then
          j=j+1
          nefluxm = nefluxm -dtest(i,j)*gradtem + vtest(i,j)*tem
        endif
        if(itport_pt(3).ne.0)then
          j=j+1
          nefluxm = nefluxm -dtest(i,j)*gradtim + vtest(i,j)*tim
        endif
        if(itport_pt(4).ne.0)then
          j=j+1
          nefluxm = nefluxm -dtest(i,j)*(nim*gradvexbm+gradnim*vexbm) 
     >     + vtest(i,j)*nim*vexbm
        endif
        if(itport_pt(5).ne.0)then
          j=j+1
          nefluxm = nefluxm -dtest(i,j)*(nim*gradvpolm+gradnim*vpolm)
     >     + vtest(i,j)*nim*vpolm
        endif
      endif     
      if(itport_pt(2).ne.0)then
        i=i+1
        j=0
        if(itport_pt(1).ne.0)then
          j=j+1
          tefluxm = tefluxm -dtest(i,j)*gradnem + vtest(i,j)*nem
        endif
        if(itport_pt(2).ne.0)then
          j=j+1
          tefluxm = tefluxm -dtest(i,j)*gradtem + vtest(i,j)*tem
        endif
        if(itport_pt(3).ne.0)then
          j=j+1
          tefluxm = tefluxm -dtest(i,j)*gradtim + vtest(i,j)*tim
        endif
        if(itport_pt(4).ne.0)then
          j=j+1
          tefluxm = tefluxm -dtest(i,j)*(nim*gradvexbm+gradnim*vexbm) 
     >    + vtest(i,j)*nim*vexbm
        endif
        if(itport_pt(5).ne.0)then
          j=j+1
          tefluxm = tefluxm -dtest(i,j)*(nim*gradvpolm+gradnim*vpolm)
     >    + vtest(i,j)*nim*vpolm
        endif
      endif     
      if(itport_pt(3).ne.0)then
        i=i+1
        j=0
        if(itport_pt(1).ne.0)then
          j=j+1
          tifluxm = tifluxm -dtest(i,j)*gradnem + vtest(i,j)*nem
        endif
        if(itport_pt(2).ne.0)then
          j=j+1
          tifluxm = tifluxm -dtest(i,j)*gradtem + vtest(i,j)*tem
        endif
        if(itport_pt(3).ne.0)then
          j=j+1
          tifluxm = tifluxm -dtest(i,j)*gradtim + vtest(i,j)*tim
        endif
        if(itport_pt(4).ne.0)then
          j=j+1
          tifluxm = tifluxm -dtest(i,j)*(nim*gradvexbm+gradnim*vexbm)
     >     + vtest(i,j)*nim*vexbm
        endif
        if(itport_pt(5).ne.0)then
          j=j+1
          tifluxm = tifluxm -dtest(i,j)*(nim*gradvpolm+gradnim*vpolm) 
     >     + vtest(i,j)*nim*vpolm
        endif
      endif     
      if(itport_pt(4).ne.0)then
        i=i+1
        j=0
        if(itport_pt(1).ne.0)then
          j=j+1
          vphifluxm = vphifluxm -dtest(i,j)*gradnem + vtest(i,j)*nem
        endif
        if(itport_pt(2).ne.0)then
          j=j+1
          vphifluxm = vphifluxm -dtest(i,j)*gradtem + vtest(i,j)*tem
        endif
        if(itport_pt(3).ne.0)then
          j=j+1
          vphifluxm = vphifluxm -dtest(i,j)*gradtim + vtest(i,j)*tim
        endif
        if(itport_pt(4).ne.0)then
          j=j+1
          vphifluxm = vphifluxm-dtest(i,j)*(nim*gradvexbm+gradnim*vexbm) 
     >    + vtest(i,j)*nim*vexbm
        endif
        if(itport_pt(5).ne.0)then
          j=j+1
          vphifluxm = vphifluxm-dtest(i,j)*(nim*gradvpolm+gradnim*vpolm)
     >    + vtest(i,j)*nim*vpolm
        endif
      endif 
      if(itport_pt(5).ne.0)then
        i=i+1
        j=0
        if(itport_pt(1).ne.0)then
          j=j+1
          vparfluxm = vparfluxm -dtest(i,j)*gradnem + vtest(i,j)*nem
        endif
        if(itport_pt(2).ne.0)then
          j=j+1
          vparfluxm = vparfluxm -dtest(i,j)*gradtem + vtest(i,j)*tem
        endif
        if(itport_pt(3).ne.0)then
          j=j+1
          vparfluxm = vparfluxm -dtest(i,j)*gradtim + vtest(i,j)*tim
        endif
        if(itport_pt(4).ne.0)then
          j=j+1
          vparfluxm = vparfluxm-dtest(i,j)*(nim*gradvexbm+gradnim*vexbm) 
     >    + vtest(i,j)*nim*vexbm
        endif
        if(itport_pt(5).ne.0)then
          j=j+1
          vparfluxm = vparfluxm-dtest(i,j)*(nim*gradvpolm+gradnim*vpolm)
     >    + vtest(i,j)*nim*vpolm
        endif
      endif   
c
      RETURN
      END  ! test_dv
