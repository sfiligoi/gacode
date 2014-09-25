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
      real*8 gradvphim
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
      gamma_p_gf = -(cv/csdam)*drhodr(jm)*(gradvexbm +gradvparm)
c
c      if(jm.eq.ngrid-1.and.itport_pt(5).eq.0)gamma_e_gf=0.0
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
       vparflux_glf = vphiflux_glf*c_per(jm)/rmajor_exp
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
