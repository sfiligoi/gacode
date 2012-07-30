       subroutine timeprofiles_dv
************************************************************************
cmnt   fills  time profiles
************************************************************************
c
      include '../inc/ptor.m'
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
c
       do j=0,mxgrid
        ne_t(j,ntime_t)=ne_m(j)
        ni_t(j,ntime_t)=ni_m(j)
        te_t(j,ntime_t)=te_m(j)
        ti_t(j,ntime_t)=ti_m(j)
        vpol_t(j,ntime_t)=vpol_m(j)
        vexb_t(j,ntime_t)=vexb_m(j)
        vphi_t(j,ntime_t)=vphi_m(j)
        doppler_shear_t(j,ntime_t)=doppler_shear_m(j)
        egamma_t(j,ntime_t)=egamma_m(j)
        anratem_t(j,ntime_t)=anrate_m(j)
        powe_t(j,ntime_t)=powe_exp(j)
        powi_t(j,ntime_t)=powi_exp(j)
       enddo
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
       return
       end  
