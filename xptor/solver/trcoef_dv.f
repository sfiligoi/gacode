ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trcoef_dv(dt)
c
c    Full dv method for the transport coefficients of VTRANS.
c    Linearises the transport equation by taking variations in both
c    the fields and field gradients. The partial derivatives of the 
c    model fluxes and sources are computed.
c
c    if iparam_pt(1) =  1 : compute fluxes, and sources over full grid
c                          transport coefficients over masked subgrid
c                    =  2 : same as 1 but do not compute convective variations
c                    = -1 : compute fluxes and sources only over full grid
c                    = -2 : compute fluxes only over masked subgrid
c
      implicit none
c
      include 'mpif.h'
      include '../inc/input.m'
      include '../inc/ptor.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/glf.m'
      include '../inc/vtrans.m'
c
      integer i,j,k,jj,sendto,nsent,k_ret,j_ret,ndone
      integer nsum
      integer m_proc,ncalls,nloops
      integer MPI_status(MPI_STATUS_SIZE)
      real*8 dne,dte,dti,dvexb,dvpol,dx,dvmin
      real*8 thetam,qm,rmajm,rminm,rhom
      real*8 tia,nia,fia,tea,nea
      real*8 tiw,niw,tew,new,nzw,testeta
      real*8 lnlam,ca,ap,am,dum,dumt,dumx,wt
      real*8 pgasa, pgasz, bgasa, bgasz
      real*8 zbrac, dt, numax, lamda, lamdamin
      real*8 ugradte,ugradti,ugradne,ugradni,ugradnz
      real*8 ugradfi,ugradfz,cnc
      real*8 fcm,r_eps,tauim,nu_plateau
      real*8 nu_braginski,vthi
      real*8 vthe,vthz
      real*8 T_ave,D_ave,x13
      real*8 mass_density,kpol1,kpol2
      real*8 last_lhs
      real*8 d_vdia,d_vneo,d_vexb,d_vpar
      real*8 gradvphim
      real*8 hm, gradhm
      real*8 unc(0:mxgrd)
      real*8 S_ei(0:mxgrd), S_pol(0:mxgrd)
      real*8 egamma_sum(0:mxgrd-1)
      real*8 gamma_p_sum(0:mxgrd-1)
      real*8 anrate_sum(0:mxgrd-1),anfreq_sum(0:mxgrd-1)
      real*8 diffgb_sum(0:mxgrd-1)
      real*8 chiegb_sum(0:mxgrd-1), chiigb_sum(0:mxgrd-1)
      real*8 chiegb_e_sum(0:mxgrd-1),cgb_sum(0:mxgrd-1)
      real*8 rhosda_sum(0:mxgrd-1),csda_sum(0:mxgrd-1)
      real*8 chiineo_sum(0:mxgrd-1),chieneo_sum(0:mxgrd-1)
      real*8 etagb_phi_sum(0:mxgrd-1)
      real*8 kpol_m_sum(0:mxgrd-1),nu_pol_m_sum(0:mxgrd-1)
      real*8 flowe_neo_sum(0:mxgrd-1)
      real*8 flowi_neo_sum(0:mxgrd-1)
      real*8 flowz_neo_sum(0:mxgrd-1)
      real*8 powe_neo_sum(0:mxgrd-1)
      real*8 powi_neo_sum(0:mxgrd-1)
      real*8 powz_neo_sum(0:mxgrd-1)
      real*8 stress_tor_i_neo_sum(0:mxgrd-1)
      real*8 stress_tor_z_neo_sum(0:mxgrd-1)
      real*8 stress_par_i_neo_sum(0:mxgrd-1)
      real*8 stress_par_z_neo_sum(0:mxgrd-1)
      real*8 flowe_adhoc_sum(0:mxgrd-1)
      real*8 flowi_adhoc_sum(0:mxgrd-1)
      real*8 flowz_adhoc_sum(0:mxgrd-1)
      real*8 powe_adhoc_sum(0:mxgrd-1)
      real*8 powi_adhoc_sum(0:mxgrd-1)
      real*8 powz_adhoc_sum(0:mxgrd-1)
      real*8 stress_tor_i_adhoc_sum(0:mxgrd-1)
      real*8 stress_tor_z_adhoc_sum(0:mxgrd-1)
      real*8 stress_par_i_adhoc_sum(0:mxgrd-1)
      real*8 stress_par_z_adhoc_sum(0:mxgrd-1)
      real*8 flowe_glf_sum(0:mxgrd-1)
      real*8 flowi_glf_sum(0:mxgrd-1)
      real*8 flowz_glf_sum(0:mxgrd-1)
      real*8 powe_glf_sum(0:mxgrd-1)
      real*8 powi_glf_sum(0:mxgrd-1)
      real*8 powz_glf_sum(0:mxgrd-1)
      real*8 stress_tor_i_glf_sum(0:mxgrd-1)
      real*8 stress_tor_z_glf_sum(0:mxgrd-1)
      real*8 stress_par_i_glf_sum(0:mxgrd-1)
      real*8 stress_par_z_glf_sum(0:mxgrd-1)
      real*8 flowe_m_sum(0:mxgrd-1)
      real*8 flowi_m_sum(0:mxgrd-1)
      real*8 flowz_m_sum(0:mxgrd-1)
      real*8 powe_m_sum(0:mxgrd-1)
      real*8 powi_m_sum(0:mxgrd-1)
      real*8 powz_m_sum(0:mxgrd-1)
      real*8 stress_tor_i_m_sum(0:mxgrd-1)
      real*8 stress_tor_z_m_sum(0:mxgrd-1)
      real*8 stress_par_i_m_sum(0:mxgrd-1)
      real*8 stress_par_z_m_sum(0:mxgrd-1)
      real*8 exch_glf_sum(0:mxgrd-1)
      real*8 xnu_m_sum(0:mxgrd-1)
      real*8 alpha_m_sum(0:mxgrd-1)
      real*8 betae_m_sum(0:mxgrd-1)
      real*8 S_ext(mxflds,mxgrd)
      real*8 glf_flux(0:10,mxflds,mxgrd),dflux(0:10,mxflds,mxgrd)
c      real*8 vdia_new_sum(3,0:mxgrd),vneo_new_sum(3,0:mxgrd)
c
      x13 = xparam_pt(13)
      cv = 1000.0
      dvmin=delt_v
c      ca = 2.D0/3.D0
      ca = 1.D0
      ap = (1.D0 + ca)/2.D0
      am = (1.D0 - ca)/4.D0
      do j=1,mxfields
      do k=1,mxgrid
        S_ext(j,k)=0.0
        do i=0,10
          glf_flux(i,j,k) = 0.D0
          dflux(i,j,k) = 0.D0
        enddo
      enddo
      enddo
      do k=0,ngrid
        unc(k)=0.D0
        S_ei(k)=0.D0
        S_pol(k)=0.D0
      enddo
      do k=0,mxgrid
       egamma_sum(k)=0.D0
       gamma_p_sum(k)=0.D0
       anrate_sum(k)=0.D0
       anfreq_sum(k)=0.D0
       diffgb_sum(k)=0.D0
       chiegb_sum(k)=0.D0
       chiigb_sum(k)=0.D0
       chiegb_e_sum(k)=0.D0
       chiineo_sum(k)=0.D0
       chieneo_sum(k)=0.0
       kpol_m_sum(k)=0.D0
       nu_pol_m_sum(k)=0.D0
       chiineogb_m(k)=0.D0
       chieneogb_m(k)=0.D0
       cgb_sum(k)=0.D0
       csda_sum(k)=0.0
       rhosda_sum(k)=0.0
       etagb_phi_sum(k)=0.D0
       flowe_neo_sum(k)=0.0
       flowi_neo_sum(k)=0.0
       flowz_neo_sum(k)=0.0
       powe_neo_sum(k)=0.0
       powi_neo_sum(k)=0.0
       powz_neo_sum(k)=0.0
       stress_tor_i_neo_sum(k)=0.0
       stress_tor_z_neo_sum(k)=0.0
       stress_par_i_neo_sum(k)=0.0
       stress_par_i_neo_sum(k)=0.0
       flowe_adhoc_sum(k)=0.0
       flowi_adhoc_sum(k)=0.0
       flowz_adhoc_sum(k)=0.0
       powe_adhoc_sum(k)=0.0
       powi_adhoc_sum(k)=0.0
       powz_adhoc_sum(k)=0.0
       stress_tor_i_adhoc_sum(k)=0.0
       stress_tor_z_adhoc_sum(k)=0.0
       stress_par_i_adhoc_sum(k)=0.0
       stress_par_z_adhoc_sum(k)=0.0
       flowe_glf_sum(k)=0.0
       flowi_glf_sum(k)=0.0
       flowz_glf_sum(k)=0.0
       powe_glf_sum(k)=0.0
       powi_glf_sum(k)=0.0
       powz_glf_sum(k)=0.0
       stress_tor_i_glf_sum(k)=0.0
       stress_tor_z_glf_sum(k)=0.0
       stress_par_i_glf_sum(k)=0.0
       stress_par_z_glf_sum(k)=0.0
       exch_glf_sum(k)=0.0
       flowe_m_sum(k)=0.0
       flowi_m_sum(k)=0.0
       flowz_m_sum(k)=0.0
       powe_m_sum(k)=0.0
       powi_m_sum(k)=0.0
       powz_m_sum(k)=0.0
       stress_tor_i_m_sum(k)=0.0
       stress_tor_z_m_sum(k)=0.0
       stress_par_i_m_sum(k)=0.0
       stress_par_z_m_sum(k)=0.0
       xnu_m_sum(k)=0.0
       alpha_m_sum(k)=0.0
       betae_m_sum(k)=0.0
       egamma_m(k)=0.D0
       anrate_m(k)=0.D0
       anfreq_m(k)=0.D0
       gamma_p_m(k)=0.D0
       nu_pol_m(k)=0.D0
       kpol_m(k)=0.D0
       csda_m(k)=0.0
       rhosda_m(k)=0.0
       cgyrobohm_m(k)=0.D0
       chiegb_m(k)=0.D0
       chiigb_m(k)=0.D0
       diffgb_m(k)=0.D0
       chiegb_etg_m(k)=0.D0
       etagb_phi_m(k)=0.D0
       etagb_par_m(k)=0.D0
       etagb_per_m(k)=0.D0
       exchgb_m(k)=0.D0
       flowe_neo(k)=0.0
       flowi_neo(k)=0.0
       flowz_neo(k)=0.0
       powe_neo(k)=0.0
       powi_neo(k)=0.0
       powz_neo(k)=0.0
       stress_tor_i_neo(k)=0.0
       stress_tor_z_neo(k)=0.0
       stress_par_i_neo(k)=0.0
       stress_par_z_neo(k)=0.0
       flowe_adhoc(k)=0.0
       flowi_adhoc(k)=0.0
       flowz_adhoc(k)=0.0
       powe_adhoc(k)=0.0
       powi_adhoc(k)=0.0
       powz_adhoc(k)=0.0
       stress_tor_i_adhoc(k)=0.0
       stress_tor_z_adhoc(k)=0.0
       stress_par_i_adhoc(k)=0.0
       stress_par_z_adhoc(k)=0.0
       flowe_glf(k)=0.0
       flowi_glf(k)=0.0
       flowz_glf(k)=0.0
       powe_glf(k)=0.0
       powi_glf(k)=0.0
       powz_glf(k)=0.0
       stress_tor_i_glf(k)=0.0
       stress_tor_z_glf(k)=0.0
       stress_par_i_glf(k)=0.0
       stress_par_z_glf(k)=0.0
       exch_glf(k)=0.0
       flowe_m(k)=0.0
       flowi_m(k)=0.0
       flowz_m(k)=0.0
       powe_m(k)=0.0
       powi_m(k)=0.0
       powz_m(k)=0.0
       stress_tor_i_m(k)=0.0
       stress_tor_z_m(k)=0.0
       stress_par_i_m(k)=0.0
       stress_par_z_m(k)=0.0
       xnu_m(k)=0.0
       alpha_m(k)=0.0
       betae_m(k)=0.0
c       do i=1,nspecies
c        vdia_new(i,k) = 0.0
c        vneo_new(i,k) = 0.0
c        vdia_new_sum(i,k) = 0.0
c        vneo_new_sum(i,k) = 0.0
c       enddo
      enddo
c
      dne = 1.D0
      dte = 1.D0
      dti = 1.D0
      dvexb = 1.D0
      dvpol = 1.D0
c      dx = arho_exp/4.D0
      dx = 1.D0
      cnc = -1.0*alpha_dia/bt_exp
c      cnc=0.0
c
c start of main dv-method loop
c
c      do k=1+i_proc,ngrid-1,n_proc
      ncalls = 2*nfields+1
      if(iparam_pt(1).eq.2)ncalls = nfields+1
      if(iparam_pt(1).le.0)ncalls = 1
      nloops = ncalls*(ngrid-1)/n_proc
c      write(*,*)"n_proc=",n_proc,"i_proc=",i_proc
c      write(*,*)"ncalls= ",ncalls," nloops = ",nloops
c      if(ncalls*(ngrid-1)/nloops.ne.n_proc)then
c       if(i_proc.eq.0)then
c        write(*,*)"warning: number of processors is not optimum"
c        write(*,*)"nearest optimum number is ",ncalls*(ngrid-1)/nloops
c       endif
c      endif
c
      do k=1,ngrid-1
        m_proc = (k-1)*ncalls
        m_proc = m_proc - n_proc*(m_proc/n_proc)
        if(iparam_pt(1).eq.-2.and.mask_r(k).eq.0)go to 21
        jm = k
        nem = (ne_m(k+1)+ne_m(k))/2.D0
c        nim = (ni_m(k+1)+ni_m(k))/2.D0
c        nzm = (nz_m(k+1)+nz_m(k))/2.D0
        tim = (ti_m(k+1)+ti_m(k))/2.D0
        tem = (te_m(k+1)+te_m(k))/2.D0
        fim = (fi_m(k+1)+fi_m(k))/2.D0
        fzm = (fz_m(k+1)+fz_m(k))/2.D0
        nim = fim*nem
        nzm = fzm*nem
        vexbm= (vexb_m(k+1)+vexb_m(k))/2.D0
        vpolm= (vpol_m(k+1)+vpol_m(k))/2.D0
        gradnem = (ne_m(k+1)-ne_m(k))/dr(k,2)
c        gradnim = (ni_m(k+1)-ni_m(k))/dr(k,2)
c        gradnzm = (nz_m(k+1)-nz_m(k))/dr(k,2)
        gradtim = (ti_m(k+1)-ti_m(k))/dr(k,2)
        gradtem = (te_m(k+1)-te_m(k))/dr(k,2)
        gradvexbm = (vexb_m(k+1)-vexb_m(k))/dr(k,2)
        gradvpolm = (vpol_m(k+1)-vpol_m(k))/dr(k,2)
        gradfim = (fi_m(k+1)-fi_m(k))/dr(k,2)
        gradfzm = (fz_m(k+1)-fz_m(k))/dr(k,2)
        gradnim = fim*gradnem + nem*gradfim
        gradnzm = fzm*gradnem + nem*gradfzm
       j_ret=0
       ipert_gf=0
       if(i_proc.eq.m_proc)then
c        write(*,*)k,"unperturbed",i_proc,j_ret
        call glf2d_dv
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
c        do i=1,nspecies
c          vdia_new(i,k) = vdia(i)
c          vneo_new(i,k) = vneo(i)
c        enddo
       endif
c
      if(iparam_pt(1).lt.1 .or. mask_r(k).eq.0)go to 21
c gradient variations
      ipert_gf=1
c      j=0
      if(itport_pt(1).ne.0)then
       j_ret = 1
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dne/dx
        gradnem = gradnem - delt_v
        gradnim = gradnim - fim*delt_v
        gradnzm = gradnzm - fzm*delt_v
c        write(*,*)"trcoef_dv",i_proc,k,j_ret
        call glf2d_dv
        gradnem = gradnem + delt_v
        gradnim = gradnim + fim*delt_v
        gradnzm = gradnzm + fzm*delt_v
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(2).ne.0)then       
       j_ret = 2
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dte/dx
        gradtem = gradtem - delt_v
c        write(*,*)k,"gradtem",i_proc,j_ret
        call glf2d_dv
        gradtem = gradtem + delt_v
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(3).ne.0)then       
       j_ret = 3
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dti/dx
        gradtim = gradtim - delt_v
c        write(*,*)k,"gradtim",i_proc,j_ret
        call glf2d_dv
        gradtim = gradtim + delt_v
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(4).ne.0)then       
       j_ret = 4
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dvexb/dx
        gradvexbm = gradvexbm - delt_v
c        write(*,*)"trcoef_dv",i_proc,k,j_ret
        call glf2d_dv
        gradvexbm = gradvexbm + delt_v
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(5).ne.0)then       
       j_ret = 5
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dvpol/dx
        gradvpolm = gradvpolm - delt_v
c        write(*,*)"trcoef_dv",i_proc,k,j_ret
        call glf2d_dv
        gradvpolm = gradvpolm + delt_v
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
c convection variations
      if(iparam_pt(1).eq.2)go to 21
c      j=0
      if(itport_pt(1).ne.0)then       
       j_ret = 6
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dne
        gradnem = gradnem*(nem+delt_v)/nem
        gradnim = gradnim*(nim+fim*delt_v)/nim
        gradnzm = gradnzm*(nzm+fzm*delt_v)/nzm
        nem = nem + delt_v
        nim = nim + fim*delt_v
        nzm = nzm + fzm*delt_v
c        write(*,*)"trcoef_dv",i_proc,k,j_ret
        call glf2d_dv
        nem = nem - delt_v
        nim = nim - fim*delt_v
        nzm = nzm - fzm*delt_v
        gradnem = gradnem*nem/(nem+delt_v)
        gradnim = gradnim*nim/(nim+fim*delt_v)
        gradnzm = gradnzm*nzm/(nzm+fzm*delt_v)
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(2).ne.0)then       
       j_ret = 7
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dte
        gradtem = gradtem*(tem+delt_v)/tem
        tem = tem + delt_v
c        write(*,*)k,"tem",i_proc,j_ret
        call glf2d_dv
        tem = tem - delt_v
        gradtem = gradtem*tem/(tem+delt_v)
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(3).ne.0)then       
       j_ret = 8
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dti
        gradtim = gradtim*(tim+delt_v)/tim
        tim = tim + delt_v
c        write(*,*)k,"tim",i_proc,j_ret
        call glf2d_dv
        tim = tim - delt_v
        gradtim = gradtim*tim/(tim+delt_v)
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(4).ne.0)then       
       j_ret = 9
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dvexb
        vexbm = vexbm + delt_v
c        write(*,*)"trcoef_dv",i_proc,k,j_ret
        call glf2d_dv
        vexbm = vexbm - delt_v
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
      if(itport_pt(5).ne.0)then       
       j_ret = 10
c       j=j+1
       m_proc = m_proc+1
       m_proc = m_proc - n_proc*(m_proc/n_proc)
       if(i_proc.eq.m_proc)then
        delt_v=dvmin*dvpol
        vpolm = vpolm + delt_v
c        write(*,*)"trcoef_dv",i_proc,k,j_ret
        call glf2d_dv
        vpolm = vpolm - delt_v
        glf_flux(j_ret,1,k) = nefluxm
        glf_flux(j_ret,2,k) = tefluxm
        glf_flux(j_ret,3,k) = tifluxm
        glf_flux(j_ret,4,k) = vphifluxm
        glf_flux(j_ret,5,k) = vparfluxm
       endif
      endif
 21    continue
      enddo
c
      call MPI_BARRIER(MPI_COMM_WORLD,i_err)
      call MPI_REDUCE(glf_flux,dflux,11*mxflds*mxgrd
     >  ,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(dflux,11*mxflds*mxgrd
     >  ,MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(egamma_m,egamma_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(egamma_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(gamma_p_m,gamma_p_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(gamma_p_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(anrate_m,anrate_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(anrate_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(anfreq_m,anfreq_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(anfreq_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(diffgb_m,diffgb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(diffgb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(chiegb_m,chiegb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(chiegb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(chiegb_etg_m,chiegb_e_sum,mxgrd,
     >  MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(chiegb_e_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(chiigb_m,chiigb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(chiigb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(chiineogb_m,chiineo_sum,mxgrd,
     >   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(chiineo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(chieneogb_m,chieneo_sum,mxgrd,
     >   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,i_err)
      call MPI_BCAST(chieneo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(cgyrobohm_m,cgb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(cgb_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(csda_m,csda_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(csda_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(rhosda_m,rhosda_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(rhosda_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(etagb_phi_m,etagb_phi_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(etagb_phi_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(kpol_m,kpol_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(kpol_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(nu_pol_m,nu_pol_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(nu_pol_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowe_neo,flowe_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowe_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowi_neo,flowi_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowi_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowz_neo,flowz_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowz_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powe_neo,powe_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powe_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powi_neo,powi_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powi_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powz_neo,powz_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powz_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_i_neo,stress_tor_i_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_i_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_z_neo,stress_tor_z_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_z_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_i_neo,stress_par_i_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_i_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_z_neo,stress_par_z_neo_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_z_neo_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowe_adhoc,flowe_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowe_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowi_adhoc,flowi_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowi_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowz_adhoc,flowz_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowz_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powe_adhoc,powe_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powe_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powi_adhoc,powi_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powi_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powz_adhoc,powz_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powz_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_i_adhoc,stress_tor_i_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_i_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_z_adhoc,stress_tor_z_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_z_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_i_adhoc,stress_par_i_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_i_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_z_adhoc,stress_par_z_adhoc_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_z_adhoc_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowe_glf,flowe_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowe_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowi_glf,flowi_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowi_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowz_glf,flowz_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowz_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powe_glf,powe_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powe_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powi_glf,powi_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powi_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powz_glf,powz_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powz_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_i_glf,stress_tor_i_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_i_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_z_glf,stress_tor_z_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_z_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_i_glf,stress_par_i_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_i_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_z_glf,stress_par_z_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_z_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(exch_glf,exch_glf_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(exch_glf_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowe_m,flowe_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowe_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowi_m,flowi_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowi_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(flowz_m,flowz_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(flowz_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powe_m,powe_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powe_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powi_m,powi_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powi_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(powz_m,powz_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(powz_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_i_m,stress_tor_i_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_i_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_tor_z_m,stress_tor_z_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_tor_z_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_i_m,stress_par_i_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_i_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(stress_par_z_m,stress_par_z_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(stress_par_z_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(xnu_m,xnu_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(xnu_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(alpha_m,alpha_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(alpha_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
      call MPI_REDUCE(betae_m,betae_m_sum,mxgrd,
     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
      call MPI_BCAST(betae_m_sum,mxgrd,MPI_DOUBLE_PRECISION
     >  ,0, MPI_COMM_WORLD, i_err)
c      call MPI_REDUCE(vdia_new,vdia_new_sum,3*mxgrd,
c     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
c      call MPI_BCAST(vdia_new_sum,mxgrd,MPI_DOUBLE_PRECISION
c     >  ,0, MPI_COMM_WORLD, i_err)
c      call MPI_REDUCE(vneo_new,vneo_new_sum,3*mxgrd,
c     > MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, i_err)
c      call MPI_BCAST(vneo_new_sum,mxgrd,MPI_DOUBLE_PRECISION
c     >  ,0, MPI_COMM_WORLD, i_err)
c
c      do i=1,nspecies
c set boundary conditions on 1/2 grid
c        vdia_new_sum(i,0)=vdia_new_sum(i,1)
c        vneo_new_sum(i,0)=vneo_new_sum(i,1)
c        vdia_new_sum(i,ngrid) = vdia_new_sum(i,ngrid-1)
c        vneo_new_sum(i,ngrid) = vneo_new_sum(i,ngrid-1)
c      enddo
c
      do k=1,ngrid-1
        egamma_m(k)=egamma_sum(k)
        gamma_p_m(k)=gamma_p_sum(k)
        anrate_m(k)=anrate_sum(k)
        anfreq_m(k)=anfreq_sum(k)
        diffgb_m(k)=diffgb_sum(k)
        chiegb_m(k)=chiegb_sum(k)
        chiigb_m(k)=chiigb_sum(k)
        chiineogb_m(k)=chiineo_sum(k)
        chieneogb_m(k)=chieneo_sum(k)
        chiegb_etg_m(k)=chiegb_e_sum(k)
        cgyrobohm_m(k)=cgb_sum(k)
        csda_m(k)=csda_sum(k)
        rhosda_m(k)=rhosda_sum(k)
        etagb_phi_m(k)=etagb_phi_sum(k)
        kpol_m(k)=kpol_m_sum(k)
        nu_pol_m(k)=nu_pol_m_sum(k)
        xnu_m(k)=xnu_m_sum(k)
        alpha_m(k)=alpha_m_sum(k)
        betae_m(k)=betae_m_sum(k)
        flowe_neo(k) = flowe_neo_sum(k)
        flowi_neo(k) = flowi_neo_sum(k)
        flowz_neo(k) = flowz_neo_sum(k)
        powe_neo(k) = powe_neo_sum(k)
        powi_neo(k) = powi_neo_sum(k)
        powz_neo(k) = powz_neo_sum(k)
        stress_tor_i_neo(k) = stress_tor_i_neo_sum(k)
        stress_tor_z_neo(k) = stress_tor_z_neo_sum(k)
        stress_par_i_neo(k) = stress_par_i_neo_sum(k)
        stress_par_z_neo(k) = stress_par_z_neo_sum(k)
        flowe_adhoc(k) = flowe_adhoc_sum(k)
        flowi_adhoc(k) = flowi_adhoc_sum(k)
        flowz_adhoc(k) = flowz_adhoc_sum(k)
        powe_adhoc(k) = powe_adhoc_sum(k)
        powi_adhoc(k) = powi_adhoc_sum(k)
        powz_adhoc(k) = powz_adhoc_sum(k)
        stress_tor_i_adhoc(k) = stress_tor_i_adhoc_sum(k)
        stress_tor_z_adhoc(k) = stress_tor_z_adhoc_sum(k)
        stress_par_i_adhoc(k) = stress_par_i_adhoc_sum(k)
        stress_par_z_adhoc(k) = stress_par_z_adhoc_sum(k)
        flowe_glf(k) = flowe_glf_sum(k)
        flowi_glf(k) = flowi_glf_sum(k)
        flowz_glf(k) = flowz_glf_sum(k)
        powe_glf(k) = powe_glf_sum(k)
        powi_glf(k) = powi_glf_sum(k)
        powz_glf(k) = powz_glf_sum(k)
        stress_tor_i_glf(k) = stress_tor_i_glf_sum(k)
        stress_tor_z_glf(k) = stress_tor_z_glf_sum(k)
        stress_par_i_glf(k) = stress_par_i_glf_sum(k)
        stress_par_z_glf(k) = stress_par_z_glf_sum(k)
        exch_glf(k) = exch_glf_sum(k)
        flowe_m(k) = flowe_m_sum(k)
        flowi_m(k) = flowi_m_sum(k)
        flowz_m(k) = flowz_m_sum(k)
        powe_m(k) = powe_m_sum(k)
        powi_m(k) = powi_m_sum(k)
        powz_m(k) = powz_m_sum(k)
        stress_tor_i_m(k) = stress_tor_i_m_sum(k)
        stress_tor_z_m(k) = stress_tor_z_m_sum(k)
        stress_par_i_m(k) = stress_par_i_m_sum(k)
        stress_par_z_m(k) = stress_par_z_m_sum(k)
c        do i=1,nspecies
c interpolate onto the full grid
c         vdia_new(i,k) = 0.5*(vdia_new_sum(i,k)+vdia_new_sum(i,k-1))
c         vneo_new(i,k) = 0.5*(vneo_new_sum(i,k)+vneo_new_sum(i,k-1))
c        enddo
      enddo
      kpol_m(0)=kpol_m(1)
      nu_pol_m(0)=nu_pol_m(1)
      xnu_m(0)=xnu_m(1)
      alpha_m(0)=0.0
      cgyrobohm_m(0) = cgyrobohm_m(1)
      csda_m(0)=csda_m(1)
      rhosda_m(0)=rhosda_m(1)
      diffgb_m(ngrid)=diffgb_m(ngrid-1)
      chiegb_m(ngrid)=chiegb_m(ngrid-1)
      chiigb_m(ngrid)=chiigb_m(ngrid-1)
      chiineogb_m(ngrid)=chiineogb_m(ngrid-1)
      chieneogb_m(ngrid)=chieneogb_m(ngrid-1)
      chiegb_etg_m(ngrid)=chiegb_etg_m(ngrid-1)
      etagb_phi_m(ngrid)=etagb_phi_m(ngrid-1)
      kpol_m(ngrid)=kpol_m(ngrid-1)
      nu_pol_m(ngrid)=nu_pol_m(ngrid-1)
      do k=ngrid,mxgrid
        cgyrobohm_m(k)=cgyrobohm_exp(k)
        csda_m(k)=csda_exp(k)
        rhosda_m(k)=rhosda_exp(k)
      enddo
c      do i=1,nspecies
c set boundary conditions on full grid
c        vdia_new(i,0)=vdia_new(i,1)
c        vneo_new(i,0)=vneo_new(i,1)
c        vdia_new(i,ngrid) = vdia_new(i,ngrid-1)
c        vneo_new(i,ngrid) = vneo_new(i,ngrid-1)
c      enddo
c
c        do k=0,jmaxm
c         write(*,155) i_proc,k,rho(k),chiineo_m(k),chiineo_sum(k)
c        enddo
c
c done with mpi part
c
c load up the matricies for vtrans
c
      lamdamin = 1.D+10
      numax = 0.D0
      do k = 1,ngrid-1
        if(iparam_pt(1).eq.-2.and.mask_r(k).eq.0)go to 11
        jm=k
        nem=(ne_m(k+1)+ne_m(k))/2.D0
c        nim=(ni_m(k+1)+ni_m(k))/2.D0
c        nzm=(nz_m(k+1)+nz_m(k))/2.D0
        fim=(fi_m(k+1)+fi_m(k))/2.D0
        fzm=(fz_m(k+1)+fz_m(k))/2.D0
        tem=(te_m(k+1)+te_m(k))/2.D0
        tim=(ti_m(k+1)+ti_m(k))/2.D0
        nim = fim*nem
        nzm = fzm*nem
        vexbm = (vexb_m(k+1)+vexb_m(k))/2.D0
        vpolm = (vpol_m(k+1)+vpol_m(k))/2.D0
c        thetam = theta_exp(k)*f_exp(k)
        thetam = theta_exp(k)
        qm=q_exp(k)
        hm =(h_m(k+1)+h_m(k))/2.0
        gradnem = (ne_m(k+1)-ne_m(k))/dr(k,2)
c        gradnim = (ni_m(k+1)-ni_m(k))/dr(k,2)
c        gradnzm = (nz_m(k+1)-nz_m(k))/dr(k,2)
        gradtem = (te_m(k+1)-te_m(k))/dr(k,2)
        gradtim = (ti_m(k+1)-ti_m(k))/dr(k,2)
        gradvexbm=(vexb_m(k+1)-vexb_m(k))/dr(k,2)
        gradvpolm=(vpol_m(k+1)-vpol_m(k))/dr(k,2)
        gradfim = (fi_m(k+1)-fi_m(k))/dr(k,2)
        gradfzm = (fz_m(k+1)-fz_m(k))/dr(k,2)
        gradnim = fim*gradnem + nem*gradfim
        gradnzm = fzm*gradnem + nzm*gradfzm
        gradhm = (h_m(k+1)-h_m(k))/dr(k,2)
c diagnostic output
c        a_unit_exp = rmin_exp(mxgrid)
        zpne_m(k) = -a_unit_exp*gradnem/nem
        zpni_m(k) = -a_unit_exp*gradnim/nim
        zpnz_m(k) = -a_unit_exp*gradnzm/nzm
        zpte_m(k) = -a_unit_exp*gradtem/tem
        zpti_m(k) = -a_unit_exp*gradtim/tim
        vthi = 9.79D3*DSQRT(2.D3*tim/amassgas_exp)/cv
        vthz = DSQRT(amassgas_exp/amassimp_exp)*vthi
        vthe = DSQRT(tem*amassgas_exp*1.8362D3/tim)*vthi
        mach_m(1,k) = (a_pol(k)*(vpol_m(k)+vneo_m(1,k))+
     >    a_tor(k)*(vexb_m(k)+vdia_m(1,k)))/vthe
        mach_m(2,k) = (a_pol(k)*(vpol_m(k)+vneo_m(2,k))+
     >    a_tor(k)*(vexb_m(k)+vdia_m(2,k)))/vthi
        mach_m(3,k) = (a_pol(k)*(vpol_m(k)+vneo_m(3,k))+
     >    a_tor(k)*(vexb_m(k)+vdia_m(3,k)))/vthz
c
c        nea = (ne_m(k)+ne_m(k-1))/2.D0
c        tia = (ti_m(k)+ti_m(k-1))/2.D0
c        fia = (fi_m(k)+fi_m(k-1))/2.D0
c        ugradti=(ti_m(k)-ti_m(k-1))/dr(k,1)
c        ugradne=(ne_m(k)-ne_m(k-1))/dr(k,1)
c        ugradfi=(fi_m(k)-fi_m(k-1))/dr(k,1)
        nea = (ne_m(k+1)+ne_m(k-1))/2.D0
        tia = (ti_m(k+1)+ti_m(k-1))/2.D0
        fia = (fi_m(k+1)+fi_m(k-1))/2.D0
        ugradti=(ti_m(k+1)-ti_m(k-1))/(dr(k,1)+dr(k,2))
        ugradne=(ne_m(k+1)-ne_m(k-1))/(dr(k,1)+dr(k,2))
        ugradfi=(fi_m(k+1)-fi_m(k-1))/(dr(k,1)+dr(k,2))
        nia = fia*nea
        ugradni = fia*ugradne+nea*ugradfi
c        nia = (ni_m(k+1)+ni_m(k-1))/2.0
c        ugradni = (ni_m(k+1)-ni_m(k-1))/(dr(k,1)+dr(k,2))
c         tiw = (ti_m(k+1)+ti_m(k)+ti_m(k-1))/3.D0
c         niw = (ni_m(k+1)+ni_m(k)+ni_m(k-1))/3.D0
c         tew = (te_m(k+1)+te_m(k)+te_m(k-1))/3.D0
c         new = (ne_m(k+1)+ne_m(k)+ne_m(k-1))/3.D0
         tew = te_m(k)
         new = ne_m(k)
         niw = ni_m(k)
         tiw = ti_m(k)
         nzw = nz_m(k)
c         niw = nia
c         tiw = tia
c mass_density = R0*(mi*ni+mz*nz)*cv , mp*10^19=1.6726D-8, mp = proton mass in kg
         mass_density = rmajor_exp*(amassgas_exp*fi_m(k)
     >                  +amassimp_exp*fz_m(k))*cv*1.6726D-8
c
c calculate electron-ion collisional energy exchange
c
c Electron-ion equilibration rates (/sec), Use NRL formula, 
c assuming Zeff is close to one, and average ion mass is amassgas_exp:
c Note that nuei is on the grid not the half grid
c Note difference due to density ne_m 1.0e13 not 1.0e14
c Note the ONETWO version is used here but not in ptorinit since have GLF here
c
        if(iexch.ne.0)then
          zbrac=(niw*zgas_exp**2 
     >    +amassgas_exp*nzw*zimp_exp**2/amassimp_exp)/new
          lnlam=24.D0-DLOG(DSQRT(1.D13*new)/(1000.D0*tew))
          lnlam=DMAX1(lnlam,1.D0)
           nuei_m(k)= 1.5D0*new*new*zbrac
     &     *3.2D-9*1.D13*lnlam
     &     /amassgas_exp/DSQRT(1000.D0*tew)**3
        endif
c
c
c  transfer fluxes
c
c       write(*,*)k,"dflux",dflux(0,2,k),dflux(2,2,k)
       nefluxm = dflux(0,1,k)
       tefluxm = dflux(0,2,k)
       tifluxm = dflux(0,3,k)
       vphifluxm = dflux(0,4,k)
       vparfluxm = dflux(0,5,k)
       neflux(k) = nefluxm
       teflux(k) = tefluxm
       tiflux(k) = tifluxm
       vphiflux(k) = vphifluxm
       vparflux(k) = vparfluxm
c
       j = 0
       if(itport_pt(1).ne.0)then
         j=j+1
         flux(j,k)=nefluxm
       endif
       if(itport_pt(2).ne.0)then
         j=j+1
         flux(j,k)=tefluxm
       endif
       if(itport_pt(3).ne.0)then
         j=j+1
         flux(j,k)=tifluxm
       endif
       if(itport_pt(4).ne.0)then
         j=j+1
         flux(j,k)=vphifluxm
       endif
       if(itport_pt(5).ne.0)then
         j=j+1
         flux(j,k)=vparfluxm
       endif
       if(mask_r(k).eq.1.and.iparam_pt(1).ge.1)then

c  find DIFF from dFLUX/dGRADT
        j = 0 
c  
        if(itport_pt(1).ne.0)then
c       transport electron density
c
         j = j+1
         delt_v=hm*dvmin*dne/dx
         nefluxm = dflux(1,1,k)
         tefluxm = dflux(1,2,k)
         tifluxm = dflux(1,3,k)
         vphifluxm = dflux(1,4,k)
         vparfluxm = dflux(1,5,k)
         i = 0
         if(itport_pt(1).ne.0)then
           i=i+1
           diff(i,j,k)=(nefluxm - neflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradnem/nem+gradhm/hm)
           vrho3(i,j,k) = (1.0/h_m(k))*1.6022D-3
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           diff(i,j,k)= (tefluxm - teflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradnem/nem+gradhm/hm)
           vrho3(i,j,k)=1.5D0*(te_m(k)/h_m(k))*1.6022D-3
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           diff(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradnem/nem+gradhm/hm)
           vrho3(i,j,k)=1.5D0*(fi_m(k)+fz_m(k))
     >                  *(ti_m(k)/h_m(k))*1.6022D-3 
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           diff(i,j,k)=(vphifluxm - vphiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradnem/nem+gradhm/hm)
           vrho3(i,j,k) = mass_density*c_tor(k)*vexb_m(k)/h_m(k)
           nu_pt(i,j,k) = -mass_density*c_tor(k)*
     >     (cnc/thetam)*tia/(dt*h_m(k))
           nu_2t(i,j,k) = -mass_density*c_tor(k)*(cnc/thetam)
     >     *(ugradti+tia*ugradfi/fia)/(dt*h_m(k))
           if(itport_pt(5).eq.0)then
             nu_2t(i,j,k) = nu_2t(i,j,k) 
     >       +mass_density*c_per(k)*(cnc/thetam)
     >        *kpol_m(k)*ugradti/(dt*h_m(k))
           else
             nu_pt(i,j,k) = nu_pt(i,j,k)*
     >      (1.0-c_per(k)*c_per(k)/(c_par(k)*c_tor(k)))
             nu_2t(i,j,k) = nu_2t(i,j,k)*
     >      (1.0-c_per(k)*c_per(k)/(c_par(k)*c_tor(k)))
           endif
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           diff(i,j,k)=(vparfluxm -vparflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradnem/nem+gradhm/hm)
         endif
       endif
c
c
       if(itport_pt(2).ne.0)then
c      transport electron energy
c
         j=j+1
         delt_v=hm*dvmin*dte/dx
         nefluxm = dflux(2,1,k)
         tefluxm = dflux(2,2,k)
         tifluxm = dflux(2,3,k)
         vphifluxm = dflux(2,4,k)
         vparfluxm = dflux(2,5,k)
         i=0 
         if(itport_pt(1).ne.0)then
           i=i+1
           diff(i,j,k)=(nefluxm - neflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtem/tem+gradhm/hm)
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           diff(i,j,k)=(tefluxm - teflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtem/tem+gradhm/hm)
           nu(i,j,k)= (nuei_m(k)/h_m(k))*1.6022D-3
           vrho3(i,j,k) = 1.5D0*(ne_m(k)/h_m(k))*1.6022D-3
           if(ABS(teflux(k)).gt.1.0D-12)then
             stiff(1,1,k) = -hm*gradtem*diff(i,j,k)/teflux(k)
           endif
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           diff(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtem/tem+gradhm/hm)
           nu(i,j,k)= - (nuei_m(k)/h_m(k))*1.6022D-3
           if(ABS(tiflux(k)).gt.1.0D-12)then
             stiff(2,1,k) = -hm*gradtem*diff(i,j,k)/tiflux(k)
           endif
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           diff(i,j,k)=(vphifluxm - vphiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtem/tem+gradhm/hm)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           diff(i,j,k)=(vparfluxm - vparflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtem/tem+gradhm/hm)
         endif
       endif
c
       if(itport_pt(3).ne.0)then
c      transport ion energy
c
         j=j+1
         delt_v=hm*dvmin*dti/dx
         nefluxm = dflux(3,1,k)
         tefluxm = dflux(3,2,k)
         tifluxm = dflux(3,3,k)
         vphifluxm = dflux(3,4,k)
         vparfluxm = dflux(3,5,k)
         i=0 
         if(itport_pt(1).ne.0)then
           i=i+1
           diff(i,j,k)=(nefluxm - neflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtim/tim+gradhm/hm)
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           diff(i,j,k)=(tefluxm - teflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtim/tim+gradhm/hm)
           nu(i,j,k)=-(nuei_m(k)/h_m(k))*1.6022D-3
           if(ABS(teflux(k)).gt.1.0D-12)then
             stiff(1,2,k) = -hm*gradtim*diff(i,j,k)/teflux(k)
           endif
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           diff(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtim/tim+gradhm/hm)
           nu(i,j,k)=(nuei_m(k)/h_m(k))*1.6022D-3
           vrho3(i,j,k)=1.5D0*(fi_m(k)+fz_m(k))
     >                  *(ne_m(k)/h_m(k))*1.6022D-3
           if(ABS(tiflux(k)).gt.1.0D-12)then
             stiff(2,2,k) = -hm*gradtim*diff(i,j,k)/tiflux(k)
           endif
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           diff(i,j,k)=(vphifluxm - vphiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtim/tim+gradhm/hm)
           nu_pt(i,j,k) = -mass_density*c_tor(k)
     >     *(cnc/thetam)*nea/(dt*h_m(k))
           nu_2t(i,j,k) = -mass_density*c_tor(k)/h_m(k)
     >     *(cnc/thetam)*(nea*ugradfi/fia+ugradne)/(dt*h_m(k))
           if(itport_pt(5).eq.0)then
             nu_pt(i,j,k) = nu_pt(i,j,k) +mass_density*c_per(k)
     >       *kpol_m(k)*(cnc/thetam)*nea/(dt*h_m(k))
           else
             nu_pt(i,j,k) = nu_pt(i,j,k)*
     >      (1.0-c_per(k)*c_per(k)/(c_par(k)*c_tor(k)))
             nu_2t(i,j,k) = nu_2t(i,j,k)*
     >       (1.0-c_per(k)*c_per(k)/(c_par(k)*c_tor(k)))
           endif
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           diff(i,j,k)=(vparfluxm - vparflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*(gradtim/tim+gradhm/hm)
           nu_p(i,j,k) = -nu_pol_m(k)*(cnc/thetam)*nea*
     >     (kpol_m(k)-(c_per(k)/c_par(k)))/h_m(k)
           nu_2(i,j,k) = -nu_pol_m(k)*(cnc/thetam)*
     >     (-c_per(k)/c_par(k))*(nea*ugradfi/fia+ugradne)/h_m(k)
         endif
       endif
c
       if(itport_pt(4).ne.0)then
c        vexb_m
c
         j=j+1
         delt_v=hm*dvmin*dvexb/dx
         nefluxm = dflux(4,1,k)
         tefluxm = dflux(4,2,k)
         tifluxm = dflux(4,3,k)
         vphifluxm = dflux(4,4,k)
         vparfluxm = dflux(4,5,k)
         i=0 
         if(itport_pt(1).ne.0)then
           i=i+1
           diff(i,j,k)=(nefluxm - neflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           diff(i,j,k)=(tefluxm - teflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           diff(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           diff(i,j,k)=(vphifluxm - vphiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
           vrho3(i,j,k) = c_tor(k)*mass_density*ne_m(k)/h_m(k)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           diff(i,j,k)=(vparfluxm - vparflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
           vrho3(i,j,k) = a_tor(k)*mass_density/h_m(k)
         endif
        endif
c
        if(itport_pt(5).ne.0)then
c         vpol_m
c
         j=j+1
         delt_v=hm*dvmin*dvpol/dx
         nefluxm = dflux(5,1,k)
         tefluxm = dflux(5,2,k)
         tifluxm = dflux(5,3,k)
         vphifluxm = dflux(5,4,k)
         vparfluxm = dflux(5,5,k)
          i=0
         if(itport_pt(1).ne.0)then
           i=i+1
           diff(i,j,k)=(nefluxm - neflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           diff(i,j,k)=(tefluxm - teflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           diff(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           diff(i,j,k)=(vphifluxm - vphiflux(k))/(delt_v)
           if(ifixeta.eq.1.and.j.eq.i)then
             if(diff(i,j,k).lt.0.0)diff(i,j,k)=ABS(eta_tor_exp(k))
           endif          
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
           vrho3(i,j,k) = c_per(k)*mass_density/h_m(k)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           diff(i,j,k)=(vparfluxm - vparflux(k))/(delt_v)
           conv3(i,j,k) = diff(i,j,k)*gradhm/hm
           nu(i,j,k) = nu_pol_m(k)/h_m(k)
           vrho3(i,j,k) = a_pol(k)*mass_density/h_m(k)
         endif
        endif
c end of diffusion terms
c
       if(mask_r(k).eq.1.and.iparam_pt(1).eq.1)then
c  calculate CONV from dFLUX/DT
        j=0
c 
        if(itport_pt(1).ne.0)then
c       transport electron density
c
         j=j+1
         delt_v=hm*dvmin*dne
         nefluxm = dflux(6,1,k)
         tefluxm = dflux(6,2,k)
         tifluxm = dflux(6,3,k)
         vphifluxm = dflux(6,4,k)
         vparfluxm = dflux(6,5,k)
         i = 0
         if(itport_pt(1).ne.0)then
           i=i+1
           conv3(i,j,k)=(nefluxm - neflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           conv3(i,j,k)=(tefluxm - teflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           conv3(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           conv3(i,j,k)=(vphifluxm - vphiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           conv3(i,j,k)=(vparfluxm - vparflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
       endif
c
       if(itport_pt(2).ne.0)then
c      transport electron energy
c
         j=j+1
         delt_v=hm*dvmin*dte
         nefluxm = dflux(7,1,k)
         tefluxm = dflux(7,2,k)
         tifluxm = dflux(7,3,k)
         vphifluxm = dflux(7,4,k)
         vparfluxm = dflux(7,5,k)
         i=0 
         if(itport_pt(1).ne.0)then
           i=i+1
           conv3(i,j,k)=(nefluxm - neflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           conv3(i,j,k)=(tefluxm - teflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           conv3(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           conv3(i,j,k)=
     >       + (vphifluxm - vphiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           conv3(i,j,k)= 
     >       + (vparfluxm - vparflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
       endif
c
       if(itport_pt(3).ne.0)then
c      transport ion energy
c
         j=j+1
         delt_v=hm*dvmin*dti
         nefluxm = dflux(8,1,k)
         tefluxm = dflux(8,2,k)
         tifluxm = dflux(8,3,k)
         vphifluxm = dflux(8,4,k)
         vparfluxm = dflux(8,5,k)
         i=0 
         if(itport_pt(1).ne.0)then
           i=i+1
           conv3(i,j,k)=(nefluxm - neflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           conv3(i,j,k)=(tefluxm - teflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           conv3(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           conv3(i,j,k)=
     >       + (vphifluxm - vphiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           conv3(i,j,k)= 
     >       + (vparfluxm - vparflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
       endif
c
       if(itport_pt(4).ne.0)then
c      transport exb velocity
c
         j=j+1
         delt_v=hm*dvmin*dvexb
         nefluxm = dflux(9,1,k)
         tefluxm = dflux(9,2,k)
         tifluxm = dflux(9,3,k)
         vphifluxm = dflux(9,4,k)
         vparfluxm = dflux(9,5,k)
         i=0 
         if(itport_pt(1).ne.0)then
           i=i+1
           conv3(i,j,k)=(nefluxm - neflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           conv3(i,j,k)=(tefluxm - teflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           conv3(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           conv3(i,j,k)= conv3(i,j,k)
     >       + (vphifluxm - vphiflux(k))/(delt_v)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           conv3(i,j,k)= conv3(i,j,k)
     >       + (vparfluxm - vparflux(k))/(delt_v)
         endif
       endif
c
       if(itport_pt(5).ne.0)then
c      transport vpol
c
         j=j+1
         delt_v=hm*dvmin*dvpol
         nefluxm = dflux(10,1,k)
         tefluxm = dflux(10,2,k)
         tifluxm = dflux(10,3,k)
         vphifluxm = dflux(10,4,k)
         vparfluxm = dflux(10,5,k)
         i=0 
         if(itport_pt(1).ne.0)then
           i=i+1
           conv3(i,j,k)=(nefluxm - neflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(2).ne.0)then
           i=i+1
           conv3(i,j,k)=(tefluxm - teflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(3).ne.0)then
           i=i+1
           conv3(i,j,k)=(tifluxm - tiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(4).ne.0)then
           i=i+1
           conv3(i,j,k)=
     >       + (vphifluxm - vphiflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
         if(itport_pt(5).ne.0)then
           i=i+1
           conv3(i,j,k)= 
     >       + (vparfluxm - vparflux(k))/(delt_v)
     >       + conv3(i,j,k)
         endif
       endif
c
c end of convection terms
        endif
c end of transport coefficients 
       endif
c
c end of main dv-method loop
c
 11   continue
c      if(iparam_pt(1).gt.0)write(*,*)"debug",k,diff(1,1,k),conv3(1,1,k)
      enddo
c
      do i=1,3
        mach_m(i,0)=mach_m(i,1)
      enddo
c
c      if(iparam_pt(1).eq.1.and.i_proc.eq.0)
c     > write(*,*)"numax = ",numax*dt,
c     >    "lamdamin = ",lamdamin/(2.D0*dr(1,1))
c reset delt_v
      delt_v=dvmin
c return for case -2 
      if(iparam_pt(1).eq.-2)return
c interpolate vrho3,nu_p,nu_2 and add time derivative terms to nu_p and nu_2
      do i=1,nfields
      do j=1,nfields
        do k=1,ngrid-1
         nu_p(i,j,k) = nu_p(i,j,k) + nu_pt(i,j,k)
         nu_2(i,j,k) = nu_2(i,j,k) + nu_pt(i,j,k)
         nu_pt(i,j,k)=0.0
         nu_2t(i,j,k)=0.0
        enddo
         vrho3(i,j,ngrid) = vrho3(i,j,ngrid-1)
         nu_p(i,j,ngrid) = nu_p(i,j,ngrid-1)
         nu_2(i,j,ngrid) = nu_2(i,j,ngrid-1)
         nu_pt(i,j,ngrid) = nu_pt(i,j,ngrid-1)
         nu_2t(i,j,ngrid) = nu_2t(i,j,ngrid-1)
         nu(i,j,ngrid) = nu(i,j,ngrid-1)
      enddo
      enddo
      if(iparam_pt(1).ge.0)then
c recompute fusion power and radiation
        if(ialpha.eq.1)then
         call palpha_dv
         if(i_proc.eq.0)write(6,*)"Pfusion = ",Pfusion
        endif
        if(irad.eq.1)then
          call prad_dv
          if(i_proc.eq.0)write(6,*)"Prad = ",Pradb_tot+Prads_tot
        endif
        if(iohm.ge.1)then
          call pohmic_dv
          if(i_proc.eq.0)write(6,*)"Pohmic = ",Pohmic_tot
        endif
      endif 
c
c compute model sources
c
      do k=2,ngrid-2
c        thetam=theta_exp(k)*f_exp(k)
        thetam=theta_exp(k)
        tim=ap*ti_m(k)+am*(ti_m(k+1)+ti_m(k-1))
        tem=ap*te_m(k)+am*(te_m(k+1)+te_m(k-1))
c        vexbm=ap*ne_m(k)*vexb_m(k)+
c     >    am*(ne_m(k+1)*vexb_m(k+1)+ne_m(k-1)*vexb_m(k-1))
        vpolm=ap*ne_m(k)*vpol_m(k)+
     >    am*(ne_m(k+1)*vpol_m(k+1)+ne_m(k-1)*vpol_m(k-1))
        S_ei(k) = nuei_m(k)*(tim-tem)*1.6022D-3
        S_pol(k) = -nu_pol_m(k)*vpolm 
c        S_pol(k) = -nu_pol_m(k)*vpolm + wp(1)*nu_pol_m(k)*unc(k)
c     >   + wp(3)*nu_pol_m(k+1)*unc(k+1)+wp(2)*nu_pol_m(k)*unc(k)
      enddo
      k=1
        thetam=theta_exp(k)
        tim=ap*ti_m(k)+am*(ti_m(k+1)+ti_m(k))
        tem=ap*te_m(k)+am*(te_m(k+1)+te_m(k))
c        vexbm=ap*ne_m(k)*vexb_m(k)+
c     >    am*(ne_m(k+1)*vexb_m(k+1)+ne_m(k)*vexb_m(k))
        vpolm=ap*ne_m(k)*vpol_m(k)+
     >    am*(ne_m(k+1)*vpol_m(k+1)+ne_m(k)*vpol_m(k))
        S_ei(k) = nuei_m(k)*(tim-tem)*1.6022D-3
        S_pol(k) = -nu_pol_m(k)*vpolm 
c        S_pol(k) = -nu_pol_m(k)*vpolm + wp(1)*nu_pol_m(k)*unc(k)
c     >   + wp(3)*nu_pol_m(k+1)*unc(k+1)
c
      k=ngrid-1
        thetam=theta_exp(k)
        tim=ap*ti_m(k)+am*(ti_m(k+1)+ti_m(k))
        tem=ap*te_m(k)+am*(te_m(k+1)+te_m(k))
c        vexbm=ap*ne_m(k)*vexb_m(k)+
c     >    am*(ne_m(k+1)*vexb_m(k+1)+ne_m(k)*vexb_m(k))
        vpolm=ap*ne_m(k)*vpol_m(k)+
     >    am*(ne_m(k+1)*vpol_m(k+1)+ne_m(k)*vpol_m(k))
        S_ei(k) = nuei_m(k)*(tim-tem)*1.6022D-3
        S_pol(k) = -nu_pol_m(k)*vpolm 
c        S_pol(k) = -nu_pol_m(k)*vpolm + wp(1)*nu_pol_m(k)*unc(k)
c     >   + wp(3)*nu_pol_m(k)*unc(k)
c
      pow_ei_cor_m(0)=0.0
      stress_par_cor_m(0)=0.0
      pow_ei_glf(0)=0.0
      do k=1,ngrid-1
        pow_ei_cor_m(k)= pow_ei_cor_m(k-1)
     >   - vprime(k,1)*dr(k,1)*S_ei(k)
        stress_par_cor_m(k) = stress_par_cor_m(k-1)
     >   - vprime(k,1)*dr(k,1)*S_pol(k)*cv*1.6726D-8
        pow_ei_glf(k) = pow_ei_glf(k-1)
     >   +vprime(k,1)*dr(k,1)*exch_glf(k)
      enddo
c
c feedback adjustments
      if(iparam_pt(1).gt.0)then
          sum_ne_m = 0.0
          sum_te_m = 0.0
          sum_ti_m = 0.0
          nsum=0
          do k=0,ngrid-1
            if(q_exp(k).gt.1.0)then
              sum_ne_m = sum_ne_m + ne_m(k)
              sum_te_m = sum_te_m + te_m(k)
              sum_ti_m = sum_ti_m + ti_m(k)
              nsum = nsum+1
            endif
          enddo
          sum_ne_m = sum_ne_m/REAL(nsum)
          sum_te_m = sum_te_m/REAL(nsum)
          sum_ti_m = sum_ti_m/REAL(nsum)
c          write(*,*)"debug",te_bc_mult,te_gain_pt,sum_te_exp
          wall_mult = wall_mult + ne_gain_pt*(sum_ne_exp-sum_ne_m)
          te_bc_mult = te_bc_mult + 
     >      te_gain_pt*(sum_te_exp-sum_te_m)/sum_te_exp
          ti_bc_mult = ti_bc_mult + 
     >      ti_gain_pt*(sum_ti_exp-sum_ti_m)/sum_ti_exp
          if(wall_mult.lt.0.0)wall_mult = 0.0
          if(te_bc_mult.lt.0.000001)te_bc_mult=0.000001
          if(ti_bc_mult.lt.0.000001)ti_bc_mult=0.000001
          if(i_proc.eq.0)then
            if(ne_gain_pt.ne.0.0)
     > write(*,*)"wall source multiplier=",wall_mult,sum_ne_exp-sum_ne_m
            if(te_gain_pt.ne.0.0)
     > write(*,*)"te bc multiplier=",te_bc_mult,sum_te_exp-sum_te_m
            if(ti_gain_pt.ne.0.0)
     > write(*,*)"ti bc multiplier=",ti_bc_mult,sum_ti_exp-sum_ti_m
          endif
      endif
c
c compute source terms (-1/vp)* d(vp*flux)/dr
c
      j=0
      if(itport_pt(1).ne.0)then
        j=j+1
        k=1
        S_ext(j,k) = wall_mult*Psour_wall(k) + smult(j)*Psour(k)
        s(j,k) = S_ext(j,k) - 
     >  (vprime(k,2)*flux(j,k))/(vprime(k,1)*dr(k,1))
c         flow_m(k) = vprime(k,2)*flux(j,k)
        INTEGRAL_RHS(j,k)= vprime(k,1)*dr(k,1)*S_ext(j,k)
        do k=2,ngrid-1
        S_ext(j,k) = wall_mult*Psour_wall(k) + smult(j)*Psour(k)
        s(j,k) = S_ext(j,k) - 
     >  (vprime(k,2)*flux(j,k)-vprime(k-1,2)*flux(j,k-1))/
     >  (vprime(k,1)*dr(k,1))
c         flow_m(k) = vprime(k,2)*flux(j,k)
         INTEGRAL_RHS(j,k) = INTEGRAL_RHS(j,k-1) + 
     >   vprime(k,1)*dr(k,1)*S_ext(j,k)
        enddo
      endif
      if(itport_pt(2).ne.0)then
        j=j+1
        k=1
        S_ext(j,k) =  smult(j)*Peaux(k)+Pe_alpha(k)+Pohpro(k)
     >     -Pradb(k)-Prads(k)
        s(j,k) = S_ext(j,k) - 
     >    (vprime(k,2)*flux(j,k))/
     >    (vprime(k,1)*dr(k,1)) 
     >    + S_ei(k)
c        powe_m(k)=vprime(k,2)*flux(j,k)
c     >   + pow_ei_cor_m(k)
        INTEGRAL_RHS(j,k) = vprime(k,1)*dr(k,1)*S_ext(j,k)
        do k=2,ngrid-1
        S_ext(j,k) =  smult(j)*Peaux(k)+Pe_alpha(k)+Pohpro(k)
     >     -Pradb(k)-Prads(k)
        s(j,k) = S_ext(j,k) - 
     >    (vprime(k,2)*flux(j,k)-vprime(k-1,2)*flux(j,k-1))/
     >    (vprime(k,1)*dr(k,1)) 
     >    + S_ei(k)
c        powe_m(k)=vprime(k,2)*flux(j,k)
c     >   + pow_ei_cor_m(k)
        INTEGRAL_RHS(j,k) = INTEGRAL_RHS(j,k-1) 
     >   + vprime(k,1)*dr(k,1)*S_ext(j,k)
        enddo
      endif
      if(itport_pt(3).ne.0)then
        j=j+1
        k=1
        S_ext(j,k) =  smult(j)*Piaux(k)+Pi_alpha(k)
        s(j,k) = S_ext(j,k) -
     >    (vprime(k,2)*flux(j,k))/
     >    (vprime(k,1)*dr(k,1)) 
     >    - S_ei(k)
c        powi_m(k)=vprime(k,2)*flux(j,k)
c     >   - pow_ei_cor_m(k)
         INTEGRAL_RHS(j,k) = vprime(k,1)*dr(k,1)*S_ext(j,k)
         do k=2,ngrid-1
        S_ext(j,k) =  smult(j)*Piaux(k)+Pi_alpha(k)
        s(j,k) = S_ext(j,k) -
     >    (vprime(k,2)*flux(j,k)-vprime(k-1,2)*flux(j,k-1))/
     >    (vprime(k,1)*dr(k,1)) 
     >    - S_ei(k)
c        powi_m(k)=vprime(k,2)*flux(j,k)
c     >   - pow_ei_cor_m(k)
         INTEGRAL_RHS(j,k) = INTEGRAL_RHS(j,k-1) +
     >   vprime(k,1)*dr(k,1)*S_ext(j,k)
cdebug        write(6,10)Piaux(k),flux(j,k),s(j,k)
c 10     format(4x,1pe12.5,4x,1pe12.5,4x,1pe12.5)
        enddo
      endif
      if(itport_pt(4).ne.0)then
c  toroidal"momentum" source in Km/s*(10^19/m^3)/s
        j=j+1
        k=1
        S_ext(j,k) = smult(j)*Mphi(k)
          s(j,k) =  S_ext(j,k) -
     >    (vprime(k,2)*flux(j,k))/
     >    (vprime(k,1)*dr(k,1))
        INTEGRAL_RHS(j,k) = vprime(k,1)*dr(k,1)*S_ext(j,k)
        do k=2,ngrid-1
          S_ext(j,k) = smult(j)*Mphi(k)
          s(j,k) = S_ext(j,k) -
     >    (vprime(k,2)*flux(j,k)-vprime(k-1,2)*flux(j,k-1))/
     >    (vprime(k,1)*dr(k,1))
          INTEGRAL_RHS(j,k) = INTEGRAL_RHS(j,k-1) +
     >    vprime(k,1)*dr(k,1)*S_ext(j,k)
        enddo
      endif
      if(itport_pt(5).ne.0)then
        j=j+1
        k=1
        S_ext(j,k) = smult(j)*Mpar(k)
        s(j,k) =  S_ext(j,k)
     > + S_pol(k)
c     > -1.D0*S_pol(k)+ wp(1)*nu_pol_m(k)*unc(k)
c     >   + nu_pol_m(k+1)*wp(3)*unc(k+1)
     >   -  vprime(k,2)*flux(j,k)/(vprime(k,1)*dr(k,1))
        INTEGRAL_RHS(j,k) = vprime(k,1)*dr(k,1)*S_ext(j,k)
        do k=2,ngrid-1
         S_ext(j,k) = smult(j)*Mpar(k)
         s(j,k) =  S_ext(j,k)
     > + S_pol(k)
c     >  -1.D0*S_pol(k) + wp(1)*nu_pol_m(k)*unc(k)
c     >   + (wp(3)*nu_pol_m(k+1)*unc(k+1)+wp(2)*nu_pol_m(k)*unc(k))
     >   - (vprime(k,2)*flux(j,k)-vprime(k-1,2)*flux(j,k-1))/
     >    (vprime(k,1)*dr(k,1))
c           if(ifixeta.eq.2.and.
c     >  (diff(nfields,nfields,k).lt.0.D0
c     >  .and.diff(nfields,nfields,k-1).lt.0.D0))then
c            if(flux(j,k)*flux(j,k-1).lt.0.D0
c     >     .and.flux(j,k)*flux(j,k+1).ge.0.D0)then
c              write(*,*)"ifixeta",k,s(j,k),s(j,k-1)
c              s(j,k) = smult(j)*Mpar(k)
c     >  -1.D0*S_pol(k) + wp(1)*nu_pol_m(k)*unc(k)
c     >     + (wp(3)*nu_pol_m(k+1)*unc(k+1)+wp(2)*nu_pol_m(k)*unc(k))
c     >     - (vprime(k+1,2)*flux(j,k+1)-vprime(k,2)*flux(j,k))/
c     >       (vprime(k,1)*dr(k,1))
c            endif
c          endif
c          stress_par_m(k) = stress_par_m(k) +stress_par_cor_m(k)
          INTEGRAL_RHS(j,k) = INTEGRAL_RHS(j,k-1) + 
     >    vprime(k,1)*dr(k,1)*S_ext(j,k)
        enddo
      endif
c
c compute the time derivative terms
c
      do i=1,nfields
        do k=2,ngrid-1
          St(i,k) = 0.D0
          do j=1,nfields
            St(i,k) = St(i,k) 
     >      + ap*vrho3(i,j,k)*(Told(j,k)-Tstart(j,k))/dt
     >      + am*vrho3(i,j,k)*(Told(j,k+1)-Tstart(j,k+1))/dt
     >      + am*vrho3(i,j,k)*(Told(j,k-1)-Tstart(j,k-1))/dt
     >      +wp(3)*(nu_2t(i,j,k+1)/2.0)
     >         *(Told(j,k+1)+Told(j,k)-Tstart(j,k+1)-Tstart(j,k))
     >      +wp(1)*(nu_2t(i,j,k)/2.0)
     >     *(Told(j,k+1)+Told(j,k-1)-Tstart(j,k+1)-Tstart(j,k-1))
     >      +wp(2)*(nu_2t(i,j,k)/2.0)
     >         *(Told(j,k)+Told(j,k-1)-Tstart(j,k)-Tstart(j,k-1))
     >      +wp(3)*(nu_pt(i,j,k+1)/dr(k,2))
     >       *(Told(j,k+1)-Told(j,k)-(Tstart(j,k+1)-Tstart(j,k)))
     >      +wp(1)*(nu_pt(i,j,k)/(dr(k,1)+dr(k,2)))
     >   *(Told(j,k+1)-Told(j,k-1)-(Tstart(j,k+1)-Tstart(j,k-1))) 
     >      +wp(2)*(nu_pt(i,j,k)/dr(k,1))
     >       *(Told(j,k)-Told(j,k-1)-(Tstart(j,k)-Tstart(j,k-1)))
          enddo
        enddo
          k=1  ! T(j,0)=T(j,1) boundary condition  
          St(i,k) = 0.D0
          do j=1,nfields
            St(i,k) = St(i,k) 
     >      + ap*vrho3(i,j,k)*(Told(j,k)-Tstart(j,k))/dt
     >      + am*vrho3(i,j,k)*(Told(j,k+1)-Tstart(j,k+1))/dt
     >      + am*vrho3(i,j,k)*(Told(j,k)-Tstart(j,k))/dt
     >      +wp(3)*(nu_2t(i,j,k+1)/2.0)
     >        *(Told(j,k+1)+Told(j,k)-Tstart(j,k+1)-Tstart(j,k))
     >      +wp(1)*(nu_2t(i,j,k)/2.0)
     >        *(Told(j,k+1)+Told(j,k)-Tstart(j,k+1)-Tstart(j,k))
     >      +wp(2)*(nu_2t(i,j,k)/2.0)
     >        *(Told(j,k)+Told(j,k)-Tstart(j,k)-Tstart(j,k))
     >      +wp(3)*(nu_pt(i,j,k+1)/dr(k,2))
     >      *(Told(j,k+1)-Told(j,k)-(Tstart(j,k+1)-Tstart(j,k)))
     >      +wp(1)*(nu_pt(i,j,k)/(dr(k,1)+dr(k,2)))
     >      *(Told(j,k+1)-Told(j,k)-(Tstart(j,k+1)-Tstart(j,k))) 
     >      +wp(2)*(nu_pt(i,j,k)/dr(k,1))
     >        *(Told(j,k)-Told(j,k)-(Tstart(j,k)-Tstart(j,k)))
          enddo
       enddo
!      dum = 0.0
!      dumt = 0.0
!      dumx = 0.0
!      do i=1,nfields
!      do k=1,ngrid-1
!         dum = dum + s(i,k)*St(i,k)
!         dumt = dumt + St(i,k)*St(i,k)
!         dumx = dumx + s(i,k)*s(i,k)
!      enddo
!      enddo
!      wt=0.0
!      if(dumt.ne.0.0)wt = dum/dumt
!      write(*,*)"check minimum",dum,dumt,wt
!      write(*,*)dumx-2.0*dum*wt+dumt*wt*wt,dumx-2.0*dum+dumt
!      if(wt.gt.100.0)then
!       if(i_proc.eq.0)write(*,*)"using newton timestep =",1.0/wt
!       do i=1,nfields
!       do k=1,ngrid-1
!         St(i,k) = wt*St(i,k)
!       enddo
!       enddo
!      endif
c transform S,St and compute integral of S
      do i=1,nfields
        last_lhs = 0.0
        do k=1,ngrid-1
          dum = S(i,k)
          S(i,k) = S(i,k) - St(i,k)
c  save the error without time derivatives 
          St(i,k)=dum 
          INTEGRAL_LHS(i,k) = last_lhs 
     >     + vprime(k,1)*dr(k,1)*(S_ext(i,k)-S(i,k))
          last_lhs = INTEGRAL_LHS(i,k)
        enddo 
      enddo
c
c      if(iparam_pt(1).gt.0)then
c       do i=1,nfields
c         dum = 0.0
c         do k=1,ngrid-1
c           dum = dum + s(i,k)*s(i,k)
c         enddo
c         dum = SQRT(dum/ngrid)
c         do k=1,ngrid-1
c           dum = ABS(s(i,k))
c           diff(i,i,k) = diff(i,i,k) + xparam_pt(12)*dr(k,1)*dum
c           diff(i,i,k) = diff(i,i,k) + 100.0
c         enddo
c       enddo
c      endif
c
c  add in numerical diffusion (Pereverzev-Corrigan)
c
      if(xparam_pt(12).ne.0.0)then
       j=0
       do i=1,5
        if(itport_pt(i).ne.0)then
         j=j+1
         do k=1,ngrid-1
          T_ave = (Told(j,k+1)+Told(j,k))/2.0
          D_ave = ABS(xparam_pt(12)*s(j,k))*dr(k,1)*dr(k,2)
          D_ave = D_ave/MAX(1.0,ABS(Told(j,k))) 
          diff(j,j,k) = diff(j,j,k) + D_ave
          conv3(j,j,k) = conv3(j,j,k) +  
     >    D_ave*((Told(j,k+1)-Told(j,k))/dr(k,2))
     >    *T_ave/MAX(1.0D-5,T_ave**2)
         enddo
        endif
       enddo 
      endif
c
 155  format('trcoef',i2,2x,i2,2x,0p1f9.6,1p6e14.6)
c
      return 
      end
