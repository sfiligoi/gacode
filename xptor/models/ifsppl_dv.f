cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ifsppl_dv
c
      Implicit None
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
      real*8 gfac,qm,zeffm,ztau_nco,zrho_nco,fac_imp_flow,art_diff
      real*8 RLT_dk, RLN_dk, RLNe_dk, q_dk, akappa_dk, shat_dk,
     &       zth_dk, xnbeam_dk, tau_dk, eps_dk, xnu_dk, rho_i_dk,
     &       v_ti_dk, rmajor_dk, g_perp_dk, RLTcrit_dk, RLTcrit2_dk,
     &       chi_0_dk, g_dk, chi_i_dk, chi_e_dk, RLTrotshear_dk,
     &       gamma_dk, rotshear_dk, rotstab
      real*8 zti
c
      cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*DSQRT(tem*1.D3)/(arho_exp*100.D0)
     >  *(1.02D2*DSQRT(tem*1.D3)/bt_exp/1.D4)**2*DSQRT(amassgas_exp)
      rhosda_m(jm)=((1.02D2*DSQRT(tem*1.D3))/bt_exp/1.D4)
     >  *DSQRT(amassgas_exp)/(arho_exp*100.D0)
      qm =(q_exp(jm+1)+q_exp(jm))/2.D0
      zeffm=(zeff_exp(jm+1)+zeff_exp(jm))/2.D0
      zpmni = -arho_exp*gradnim/nim
      zpmte=-arho_exp*gradtem/tem
      zpmti=-arho_exp*gradtim/tim
c      zti=-(dlog(ti_m(jm-1))-dlog(ti_m(jm)))/(rho(jm-1)-rho(jm))
      zpni_m(jm)=zpmni
      zpte_m(jm)=zpmte
      zpti_m(jm)=zpmti
      csdam=9.79D5*DSQRT(tem*1.D3)/
     &      (arho_exp*100.D0)/DSQRT(amassgas_exp)
      csda_m(jm)=csdam
c      write(*,*) 'zpmti = ',zpmti,zti
c
c...  IFS/PPPL (Dorland Kotchenreuther) model
c     Two versions are implemented - The 1994 original version (ip_chi1)
c     and the 1996 Varenna version (ip_chi2) with both velocity and
c     diamagnetic shear stabilization plus an additional elongation factor.
c
c       write(*,*) 'iv_dk = ',iv_dk
       RLT_dk=rmajor_exp/arho_exp*zpmti*drhodr(jm)
       zpmni=zpmne
       RLN_dk=rmajor_exp/arho_exp*zpmni*drhodr(jm)
       RLNe_dk=RLN_dk
       q_dk=q_exp(jm)
       akappa_dk=elong_exp(jm)
       shat_dk=abs(shat_exp(jm))*drhodrrrho(jm)
       zth_dk=zeff_exp(jm)
       xnbeam_dk=nfast_exp(jm)/ne_exp(jm)
       tau_dk=tim/tem
       eps_dk=rmin_exp(jm)/rmajor_exp
       xnu_dk=2.5D-7*nem*1.D13/((tem*1.D3)**1.5D0*(tim*1.D3)**0.5D0)
     >  *rmajor_exp
       rho_i_dk=rhosda_m(jm)*arho_exp*dsqrt(tim/tem)
       v_ti_dk=csda_m(jm)*arho_exp*dsqrt(tim/tem)
       rmajor_dk=rmajor_exp
       g_perp_dk=0.D0
c
c..jek Varenna 1996 version of IFS model
c      set shearing rate g_perp zero and compute
c      rotstab outside of ip_chi routine
c
       if(iv_dk.eq.2) then
         call ip_chi2(itest_ntcc,itest_dk,
     &     RLT_dk,RLN_dk,RLNe_dk,q_dk,akappa_dk,
     &     shat_dk,zth_dk,
     &     xnbeam_dk,tau_dk,eps_dk,xnu_dk,g_perp_dk,
     &     rmajor_dk,rho_i_dk,v_ti_dk,RLTcrit_dk,
     &     RLTcrit2_dk,chi_0_dk,g_dk,
     &     gamma_dk,chi_i_dk,chi_e_dk)
       endif
       gamma_dk=gamma_dk/csda_m(jm)
       if (gamma_dk.le.0.) gamma_dk=1.D-6
c
c... leave out ExB shear for now
c
       rotstab=1.D0
       if ( irotstab.eq.0 ) then
         rotshear_dk=alpha_e_dk*egamma_exp(jm)+
     &              alpha_star_dk*vstarp_exp(jm)
         rotstab=(1.D0-abs(rotshear_dk/
     &          (gamma_dk*csda_m(jm)/csda_exp(jm))))
       elseif (irotstab.gt.0) then
         rotshear_dk=alpha_e_dk*egamma_m(jm)+
     &              alpha_star_dk*vstarp_m(jm)
         rotstab=(1.D0-abs(rotshear_dk/gamma_dk))
       else
         rotshear_dk=0.D0
       endif
c
       if(igeo_m.ge.1) gfac=geofac(jm)
c
c       write(*,*) 'chii_dk = ',chi_i_dk
c       write(*,*) 'chie_dk = ',chi_e_dk
c       write(*,*) 'rotstab = ',rotstab
c       write(*,*) 'irotstab = ',irotstab
c       write(*,*) 'rotshear_dk = ',rotshear_dk
c       write(*,*) 'gamma_dk = ',gamma_dk
c       write(*,*) 'alpha_e,alpha_star  =',alpha_e_dk,alpha_star_dk
c       write(*,*) 'egamma_exp = ',egamma_exp(jm)
c       write(*,*) 'egamma_m = ',egamma_m(jm)
c       write(*,*) 'vstarp_exp = ',vstarp_exp(jm)
c       write(*,*) 'vstarp_m = ',vstarp_m(jm)
c
       chietem=cmodel_e*rotstab*chi_e_dk*gfac*cgyrobohm_m(jm)
       chiitim=cmodel_i*rotstab*chi_i_dk*gfac*cgyrobohm_m(jm)
       diffnem=0.D0
c
       zpti_dk(jm)=
     >    RLT_dk/rmaj_exp(jm)*arho_exp/drhodr(jm)
       zptim_dk(jm)=
     >    RLTcrit_dk/rmaj_exp(jm)*arho_exp/drhodr(jm)
       zptim2_dk(jm)=
     >    RLTcrit2_dk/rmaj_exp(jm)*arho_exp/drhodr(jm)
c
       chi0_dk(jm)=chi_0_dk
       g0_dk(jm)=g_dk
c
       chie_m(jm)=chietem
       chii_m(jm)=chiitim
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
c... compute fluxes, power flows
c
       tefluxm = (1.6022D-3)*tem*nem*zpmte*chietem/arho_exp
       tifluxm = (1.6022D-3)*tim*nim*zpmti*chiitim/arho_exp
c      powem = kevdsecpmw*tem*nem*1.D19/arho_exp*gradrhosq_exp(jm)*
c     &        sfactor(jm)*(chietem*zpmte)
c     &        +xconv*1.5D0*tem*flow_exp(jm)
c      powim = fac_imp_flow*
c     &        kevdsecpmw*tim*nim*1.D19/arho_exp*gradrhosq_exp(jm)*
c     &        sfactor(jm)*(chiitim*zpmti)
c     &        +xconv*1.5D0*tim*flow_exp(jm)
c      powem = 1.6022D-3*vprime(jm,2)*tefluxm + pow_ei_cor_m(jm)
c      powim = 1.6022D-3*vprime(jm,2)*tifluxm - pow_ei_cor_m(jm)
c
c      write(*,50) jm, rho(jm), powi_exp(jm), powim,
c     &            powe_exp(jm), powem
c      
 50   format(2x,i2,2x,0p1f10.6,1p6e14.6)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
