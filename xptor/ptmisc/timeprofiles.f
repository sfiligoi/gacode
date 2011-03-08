       subroutine timeprofiles(ktime)
************************************************************************
c 22-mar-00 jek added provision for artificial Guassian chie at a
c           particular time-step 
cmnt   fills  time profiles
************************************************************************
c
      include '../inc/ptor.m'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
c
       integer ktime
       real*8 chie_art(0:jmaxm)
c
       do j=0,jmaxm
        te_t(j,ktime)=te_m(j)
        ti_t(j,ktime)=ti_m(j)
        ne_t(j,ktime)=ne_m(j)
        ni_t(j,ktime)=ni_m(j)
        nz_t(j,ktime)=nz_m(j)
        vphi_t(j,ktime)=vphi_m(j)*corot
        vexb_t(j,ktime)=vexb_m(j)
        vpol_t(j,ktime)=vpol_m(j)
        vetor_t(j,ktime)=vetor_m(j)
        vepol_t(j,ktime)=vepol_m(j)
        vstar_t(j,ktime)=vstar_m(j)
        egamma_t(j,ktime)=egamma_m(j)
        egamma_vphi_t(j,ktime)=egamma_vphi(j)
        egamma_vpol_t(j,ktime)=egamma_vpol(j)
        egamma_vstar_t(j,ktime)=egamma_vstar(j)
        anratem_t(j,ktime)=anrate_m(j)
        chie_t(j,ktime)=chiegb_m(j)*cgyrobohm_m(j)
        chii_t(j,ktime)=chiigb_m(j)*cgyrobohm_m(j)
c        write(*,100) j, rho(j), chii_m(j), chiigb_m(j), chii_t(j,ktime)
        etaphi_t(j,ktime)=etagb_phi_m(j)*cgyrobohm_m(j)
        gammap_t(j,ktime)=gamma_p_m(j)
        chie_art(j)=0.D0
        if (ichiep.gt.0 .and. ktime.eq.ktime_chiep) then
          chie_art(j)=chiep_c*dexp(-chiep_a*(rho(j)-chiep_b)**2)
          chiegb_m(j)=chiegb_m(j) + chie_art(j)/cgyrobohm_m(j)
        endif 
       enddo
       if(xon_pt.lt.1.) then
           te_t(jmaxm,ktime)=1.D-10
           ti_t(jmaxm,ktime)=1.D-10
       endif
c
       if (time_series.ne.0) then
         do j=0,jmaxm
          te_exp_t(j,ktime)=te_exp(j)
          ti_exp_t(j,ktime)=ti_exp(j)
          vphi_exp_t(j,ktime)=vphi_exp(j)
         enddo
         powe_exp_t(ktime)=powe_exp(jmaxm)
         powi_exp_t(ktime)=powi_exp(jmaxm)
         powewdot_exp_t(ktime)=powe_wdot_exp(jmaxm)
         powiwdot_exp_t(ktime)=powi_wdot_exp(jmaxm)
       endif
c
 100  format(i2,2x,0p1f5.2,1p6e14.5)
       return
       end  
