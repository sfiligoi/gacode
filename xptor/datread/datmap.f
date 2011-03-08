c@datmap.f
c jek 12-feb-09
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c... Maps data from "_d" variables to "_exp" variables
c    Note: some _exp variables, jmaxm in input.m (e.g. arho_exp)
c          glf.m needed for i_proc
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  
      subroutine datmap(ialpha)
c
      implicit none
c
      include '../inc/data.m'
      include '../inc/tport.m'
      include '../inc/input.m'
      include '../inc/glf.m'
      include '../inc/model.m'
c
      integer j, jj, ialpha, ipure
      real*8 drm, dvoldr_p, dvoldr_m, shift_shaf, dummy1, aomega
      real*8 dum(nj)
c
      pi_m=atan2(0.0D0,-1.0D0)
      ipure=0
c
      if(nscale.eq.0) nscale=1.D0
c
c... 1D variable transfers
c
      arho_exp = r_d(jmaxm+1)
      rmajor_exp = rmajor_d
      bt_exp = dabs(btor_d)
      sign_bt_exp = btor_d/bt_exp
      bp0_exp = dabs(bp0_d(1))
      elonga_exp=elongx_d(jmaxm+1)
      deltaa_exp=deltax_d(jmaxm+1)
c
c... 2d profile transfers
c
      do j=1,nj_d
        te_exp(j-1)=te_d(j)
        ti_exp(j-1)=ti_d(j)
        ne_exp(j-1)=nscale*1.D-19*ene_d(j)
        zeff_exp(j-1)=zeff_d(j)
        if(ipure.eq.1) zeff_exp(j-1)=1.D0
c... only 2 max primary ions for now
        if(idata.eq.1 .and. nprim_d.eq.2) then
          ni_exp(j-1)=nscale*1.D-19*(en_d(j,1)+en_d(j,2))  ! D, T
        else
          ni_exp(j-1)=nscale*1.D-19*en_d(j,1)
        endif
c
        if(ipure.eq.1) ni_exp(j-1)=ne_exp(j-1)
c... only 2 max impurity ions
        if(idata.eq.1 .and. nprim_d.eq.2) then  ! impurity is assumed no. 3 ion
          if(nimp_d.eq.1) then
            nz_exp(j-1)=nscale*1.D-19*en_d(j,3)
          else
            nz_exp(j-1)=nscale*1.D-19*(en_d(j,3)+en_d(j,4))
          endif
        else
          nz_exp(j-1)=nscale*1.D-19*en_d(j,2)
        endif
        nitot_exp(j-1)=nscale*ni_exp(j-1)+nz_exp(j-1)
        if(ipure.eq.1) nitot_exp(j-1)=ne_exp(j-1)
        nfst_exp(j-1)=nscale*1.D-19*enbeam_d(j)
        ptot_exp(j-1)=nscale*ptot_d(j)*1.0D-19
        pfast_exp(j-1)=nscale*1.D-19*pfast_d(j)
        torque_exp(j-1)=torque_d(j)
        psir_exp(j-1)=psir_d(j)
        q_exp(j-1)=q_d(j)
        rmin_exp(j-1)=rminavnpsi_d(j)
        rmaj_exp(j-1)=rmajavnpsi_d(j)
        angrot_exp(j-1)=angrot_d(j)
        f_exp(j-1)=fcap_d(j)
        g_exp(j-1)=gcap_d(j)
        h_exp(j-1)=hcap_d(j)
        gradrhosq_exp(j-1)=grho2npsi_d(j)
        gradrho_exp(j-1)=grho1npsi_d(j)
        volume_exp(j-1)=psivolp_d(j)
        elong_exp(j-1)=elongx_d(j)
        delta_exp(j-1)=deltax_d(j)
        sfarea_exp(j-1)=sfareanpsi_d(j)
        cxarea_exp(j-1)=cxareanpsi_d(j)
        curden_exp(j-1) = curden_d(j)
        curboot_exp(j-1) = curboot_d(j)
        curbeam_exp(j-1) = curbeam_d(j)
        currf_exp(j-1) = currf_d(j)
        curohm_exp(j-1) = curohm_d(j)
        qbeame_exp(j-1) = qbeame_d(j)
        qbeami_exp(j-1) = qbeami_d(j)
        qrfe_exp(j-1) = qrfe_d(j)
        qrfi_exp(j-1) = qrfi_d(j)
        qohm_exp(j-1) = qohm_d(j)
        qrad_exp(j-1) = qrad_d(j)
        qione_exp(j-1) = qione_d(j)
        qioni_exp(j-1) = qioni_d(j)
        sbeame_exp(j-1)=sbion_d(j)
        sbeam_exp(j-1)=sbeam_d(j)
c
        do jj=1,nprim_d
          sion_exp(j-1,jj)=sion_d(j,jj)
          srecom_exp(j-1,jj)=srecom_d(j,jj)
          scx_exp(j-1,jj)=scx_d(j,jj)
          sbcx_exp(j-1,jj)=sbcx_d(j,jj)
          s_exp(j-1,jj)=s_d(j,jj)
          sdot_exp(j-1,jj)=dudtsv_d(j,jj)
        enddo
        do jj=1,nneu_d
          enn_exp(j-1,jj)=enn_d(j,jj)
          ennw_exp(j-1,jj)=ennw_d(j,jj)
          ennv_exp(j-1,jj)=ennv_d(j,jj)
          volsn_exp(j-1,jj)=volsn_d(j,jj)
        enddo
c
        if(ncl_flag.eq.1)then
          xb2_exp(j-1)=xb2_d(j)
          xbm2_exp(j-1)=xbm2_d(j)
          xngrth_exp(j-1)=xngrth_d(j)
          xgrbm2_exp(j-1)=xgrbm2_d(j)
          fm1_exp(j-1)=fm1_d(j)
          fm2_exp(j-1)=fm2_d(j)
          fm3_exp(j-1)=fm3_d(j)
          fhat_exp(j-1)=fhat_d(j)
        endif
c... some ufile data
        nm1_exp(j-1)=nscale*1.D-19*en_nm1_d(j)
        nm2_exp(j-1)=nscale*1.D-19*en_nm2_d(j)
        nm3_exp(j-1)=nscale*1.D-19*en_nm3_d(j)
      enddo
c smooth rmin gradient
      do j=1,nj_d-1
        dum(j+1) = rmin_exp(j)-rmin_exp(j-1)
      enddo
      dum(1) = dum(2)
      if(igyro.ne.1) call average7_1d(dum,nj_d)
      rmin_exp(0)=0.0
      do j=1,nj_d-1
        rmin_exp(j) = rmin_exp(j-1)+dum(j+1)
      enddo
c smooth rmaj gradient
      do j=1,nj_d-1
        dum(j+1) = rmaj_exp(j)-rmaj_exp(j-1)
      enddo
      dum(1) = dum(2)
      if(igyro.ne.1) call average7_1d(dum,nj_d)
      do j=1,nj_d-1
        rmaj_exp(j) = rmaj_exp(j-1)+dum(j+1)
      enddo
! smooth delta gradient
      do j=1,nj_d-1
        dum(j+1) = delta_exp(j)-delta_exp(j-1)
      enddo
      dum(1) = dum(2)
      if(igyro.ne.1) call average7_1d(dum,nj_d)
      do j=1,nj_d-1
        delta_exp(j) = delta_exp(j-1)+dum(j+1)
      enddo
! smooth elong gradient
      do j=1,nj_d-1
        dum(j+1) = elong_exp(j)-elong_exp(j-1)
      enddo
      dum(1) = dum(2)
      if(igyro.ne.1) call average7_1d(dum,nj_d)
      do j=1,nj_d-1
        elong_exp(j) = elong_exp(j-1)+dum(j+1)
      enddo
c
c... 2d energy and particle flows
c    central values needed up front
      vol_exp(0)=0.D0
      powe_beam_exp(0)=0.D0
      powe_rf_exp(0)=0.D0
      powe_lh_exp(0)=0.D0
      powe_oh_exp(0)=0.D0
      powe_rad_exp(0)=0.D0
      powe_ion_exp(0)=0.D0
      powe_wdot_exp(0)=0.D0
      powe_fus_exp(0)=0.D0
      powi_beam_exp(0)=0.D0
      powi_rf_exp(0)=0.D0
      powi_ion_exp(0)=0.D0
      powi_cx_exp(0)=0.D0
      powi_wdot_exp(0)=0.D0
      powi_fus_exp(0)=0.D0
      pow_ei_exp(0)=0.D0
      powe_exp(0)=0.D0
      powi_exp(0)=0.D0
      flow_wall_exp(0)=0.D0
      flow_recom_exp(0)=0.D0
      flow_beam_exp(0)=0.D0
      flow_sdot_exp(0)=0.D0
      flow_exch_exp(0)=0.D0
      flow_exp(0)=0.D0  
      if(ascale.eq.0.) ascale=1.D0
      if(bscale.eq.0.) bscale=1.D0
      bt_exp = abs(btor_d)*bscale
      if(pscale.eq.0) pscale=1.D0
      if(pbescale.eq.0) pbescale=1.D0
      if(pbiscale.eq.0) pbiscale=1.D0
      if(snbscale.eq.0) snbscale=1.D0
      if(nfscale.eq.0) nfscale=1.D0
      if(prfscale.eq.0) prfscale=1.D0
      if(prfescale.eq.0) prfescale=1.D0
      if(prfiscale.eq.0) prfiscale=1.D0
      if(wallneut.eq.0) wallneut=1.D0
      if(wallneutp.eq.0) wallneutp=1.D0
c
      do j=2,nj_d    
        drm=r_d(j)-r_d(j-1)
        dvoldr_p=2.D0*pi_m*r_d(j)*2.D0*pi_m*rmajor_d*hcap_d(j)
c        write(*,'(i2,2x,0pf10.5,1p2e12.4)')j,rho_d(j),r_d(j),hcap_d(j)
        dvoldr_m=2.D0*pi_m*r_d(j-1)*2.D0*pi_m*rmajor_d*hcap_d(j-1)
        vol_exp(j-1)=vol_exp(j-2)+
     >       0.5D0*(dvoldr_p+dvoldr_m)*drm
c power in MW
        powe_beam_exp(j-1)=powe_beam_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qbeame_d(j)+dvoldr_m*qbeame_d(j-1))*drm
        powe_rf_exp(j-1)=powe_rf_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qrfe_d(j)+dvoldr_m*qrfe_d(j-1))*drm
        powe_lh_exp(j-1)=powe_lh_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qlhe_d(j)+dvoldr_m*qlhe_d(j-1))*drm
        powe_oh_exp(j-1)=powe_oh_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qohm_d(j)+dvoldr_m*qohm_d(j-1))*drm
        powe_rad_exp(j-1)=powe_rad_exp(j-2)+ 
     >   1.D-6*0.5D0*(dvoldr_p*qrad_d(j)+dvoldr_m*qrad_d(j-1))*drm
        powe_ion_exp(j-1)=powe_ion_exp(j-2)+wallneutp*
     >   1.D-6*0.5D0*(dvoldr_p*qione_d(j)+ dvoldr_m*qione_d(j-1))*drm
        powe_wdot_exp(j-1)=powe_wdot_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*dpedtc_d(j)+dvoldr_m*dpedtc_d(j-1))*drm
        powe_fus_exp(j-1)=powe_fus_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qfuse_d(j)+dvoldr_m*qfuse_d(j-1))*drm
        powi_beam_exp(j-1)=powi_beam_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qbeami_d(j)+dvoldr_m*qbeami_d(j-1))*drm
        powi_rf_exp(j-1)=powi_rf_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qrfi_d(j)+dvoldr_m*qrfi_d(j-1))*drm
        powi_ion_exp(j-1)=powi_ion_exp(j-2)+wallneutp*
     >   1.D-6*0.5D0*(dvoldr_p*qioni_d(j)+dvoldr_m*qioni_d(j-1))*drm
        powi_cx_exp(j-1)=powi_cx_exp(j-2)+wallneutp*
     >   1.D-6*0.5D0*(dvoldr_p*qcx_d(j)+dvoldr_m*qcx_d(j-1))*drm
        powi_wdot_exp(j-1)=powi_wdot_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*dpidtc_d(j)+dvoldr_m*dpidtc_d(j-1))*drm
        powi_fus_exp(j-1)=powi_fus_exp(j-2)+
     >   1.D-6*0.5D0*(dvoldr_p*qfusi_d(j)+dvoldr_m*qfusi_d(j-1))*drm
        pow_ei_exp(j-1)=pow_ei_exp(j-2)-
     >   1.D-6*0.5D0*(dvoldr_p*qdelt_d(j)+dvoldr_m*qdelt_d(j-1))*drm
c flows in KA
        flow_wall_exp(j-1)=flow_wall_exp(j-2)+wallneut*
     >       1.6022D-21*0.5D0*(dvoldr_p*(sion_d(j,1))+
     >       dvoldr_m*(sion_d(j-1,1)))*drm
        flow_recom_exp(j-1)=flow_recom_exp(j-2)+
     >       1.6022D-21*0.5D0*(dvoldr_p*(srecom_d(j,1))+
     >       dvoldr_m*(srecom_d(j-1,1)))*drm 
        flow_beam_exp(j-1)=flow_beam_exp(j-2)+1.6022D-21*0.5D0*
     >       (dvoldr_p*(snbscale*sbeam_d(j)+sbcx_d(j,1))+
     >       dvoldr_m*(snbscale*sbeam_d(j-1)+sbcx_d(j-1,1)))*drm
        flow_sdot_exp(j-1)=flow_sdot_exp(j-2)+
     >       1.6022D-21*0.5D0*(dvoldr_p*dudtsv_d(j,1)+
     >       dvoldr_m*dudtsv_d(j-1,1))*drm
      enddo
c
c... total flows
c    Note: switches (xoh_exp,xwdot,xfus_exp) in input.m
c
      if(idata.eq.1) then
        xrad_exp=-xrad_exp
        xion_exp=-xion_exp
      endif
c
c      write(*,*) 'xfus_exp = ',xfus_exp
      do j=2,nj_d
        powe_exp(j-1)=pbescale*powe_beam_exp(j-1)+
     >   prfscale*prfescale*powe_rf_exp(j-1)+powe_lh_exp(j-1)+
     >   (1.D0-xoh_exp)*powe_oh_exp(j-1)
     >   -xrad_exp*powe_rad_exp(j-1)-xion_exp*powe_ion_exp(j-1)-
     >   (1.D0-xwdot)*powe_wdot_exp(j-1)
     >   -pow_ei_exp(j-1)
     >   +xfus_exp*powe_fus_exp(j-1)
        powi_exp(j-1)=pbiscale*powi_beam_exp(j-1)+
     >   prfscale*prfiscale*powi_rf_exp(j-1)
     >   -xion_exp*powi_ion_exp(j-1)+powi_cx_exp(j-1)-
     >   (1.D0-xwdot)*powi_wdot_exp(j-1)
     >   +pow_ei_exp(j-1)
     >   +xfus_exp*powi_fus_exp(j-1)
        flow_exp(j-1)=flow_wall_exp(j-1)+flow_recom_exp(j-1)+
     >                nfscale*flow_beam_exp(j-1)-
     >                (1.D0-xsdot)*flow_sdot_exp(j-1)
      enddo
c
c... Exchange flow
c    flow_exch_exp(j) is the expansion cooling in the ion channel integrated
c    to a power flow (MW); negative ie cooling for positive plasma flow.
c    It estimates the possible size of "anomalous" energy exchange. It should
c    be compared with pow_ei_exp(j).For electron directed waves we expect 
c    electron channel cooling and ion channel heating by this amount.      
c    For ion directed waves we expect electron channel heating and ion 
c    channel cooling by this amount
c
       flow_exch_exp(0)=0.D0
       do j=1,nj_d-1
        if(j.eq.1)then
          drm = r_d(2)-r_d(1)
        else
          drm=r_d(j)-r_d(j-1)
        endif
        flow_exch_exp(j)=flow_exch_exp(j-1)+
     >  (kevdsecpmw/1.0D-19)*
     >       drm*flow_exp(j)*ti_exp(j)*(dlog(ti_exp(j)*ni_exp(j))
     >       -dlog(ti_exp(j-1)*ni_exp(j-1)))/drm
       enddo
c
c... rescale total power flows if desired
c
       do j=2,nj_d
        powi_exp(j-1)=powi_exp(j-1)*pscale
        powe_exp(j-1)=powe_exp(j-1)*pscale
        pow_ei_exp(j-1)=pow_ei_exp(j-1)*pscale
       enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... Neoclassical transport transfers
c
      if (use_xneo_m.eq.-1) then
        do j=1,nj_d
          chiineo_exp(j-1)=xkineo_d(j)/gradrhosq_exp(j-1)
        enddo
      endif
      if (use_xneo_m.le.1 .and. use_xneo_m.ge.0) then
        do j=1,nj_d
          chiineo_exp(j-1)=xkineo_d(j)/ni_exp(j-1)/1.D19/
     >                     (bscale*ascale**(.5D0))
     >                     /gradrhosq_exp(j-1)
          chieneo_exp(j-1)=chiineo_exp(j-1)/42.D0/
     >                     dsqrt(amassgas_exp)
        enddo
      endif
      if(ncl_flag.eq.0) xb2_exp(0)=xb2_exp(1)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c  surface factor
c
      do j=0,jmaxm
         sfactor(j)=2.D0*pi_m*arho_exp*rho(j)*
     >              h_exp(j)*2.D0*pi_m*rmajor_exp
c       write(*,'(i2,2x,0p6f10.5)') j, rho(j), sfactor(j)
      enddo
c
c...noncircular geometric factors
c
      do j=1,jmaxm-1
        drhodr(j)=(rho(j+1)-rho(j))*arho_exp/
     >   (rmin_exp(j+1)-rmin_exp(j))
        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/arho_exp/rho(j)
        shift_shaf=(rmaj_exp(j)-rmaj_exp(j-1))/
     >   (rmin_exp(j)-rmin_exp(j-1))
        geofac(j)=gradrho_exp(j)*(rho(j)-rho(j-1))*arho_exp
     >   /(rmin_exp(j)-rmin_exp(j-1))/gradrhosq_exp(j)
        georotrate(j)=elong_exp(j)/(1.D0+shift_shaf)**2
        geoalpha(j)=1.D0/rmajor_exp*(vol_exp(j)-vol_exp(j-1))
     >              /(rho(j)-rho(j-1))/arho_exp/
     >              (2.D0*pi_m*rho(j)*arho_exp)**2*
     >              dsqrt(vol_exp(j)/(2.D0*pi_m**2*rmajor_exp))
        bteff_exp(j)=bt_exp*rho(j)*arho_exp/rmin_exp(j)*drhodr(j)
      enddo
c
c... magnetic shear
c
      do j=1,jmaxm-1
         shat_exp(j)=(rho(j+1)+rho(j))/(q_exp(j+1)+q_exp(j))*
     >               (q_exp(j+1)-q_exp(j))/(rho(j+1)-rho(j))
         if (DABS(shat_exp(j)).LT.1.D-10) shat_exp(j)=1.D-10
      enddo
c
c... line-average electron density (10^19 m^-3)
c
      call trapl(r_d*100.D0,ene_d*1.D-19,nj_d,nebar_exp)
      nebar_exp=nebar_exp/(r_d(nj_d)*100.D0)
c
c... if grad-rho is missing
c
      if (grho1npsi_d(2).eq.0) then
         if (i_proc.eq.0) write(*,80)
         do j=1,nj_d
           gradrhosq_exp(j-1)=(1.D0+elonga_exp**2)/2.D0/elonga_exp
           gradrho_exp(j-1)=dsqrt(gradrhosq_exp(j-1))
         enddo
      endif
c
      jstinv=1
      do j=2,jmaxm
        if (q_exp(j).gt.1..and.q_exp(j-1).le.1.) jstinv=j
      enddo
      jstmix=dsqrt(2.D0)*jstinv
c
c... Experimental Er profile
c
      do j=1,jmaxm
        er_exp(j-1)=er_d(j)
      enddo
c
c      do j=0,jmaxm
c        write(*,'i2,2x,0p6f10.5') j, rho(j), er_exp(j)
c      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... read in experimental ExB shear profile
c
      if(iexb.eq.1) then
        open(unit=5,status='unknown',access='sequential',
     >             file='omega.dat')
        do j=0,jmaxm
          csda_exp(j)=9.79D5*dsqrt(te_exp(j)*1.D3)/
     >        (arho_exp*100.D0)/dsqrt(amassgas_exp)
          read(5,*) dummy1,aomega
          megamma_exp(j)=1.D3*aomega/csda_exp(j)
        enddo
        close(5)
      endif

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... set central,edge values (overwrites)
c
      shat_exp(0)=0.D0
      rmin_exp(0)=1.D-6
      drhodr(0)=drhodr(1)
      drhodrrrho(0)=drhodrrrho(1)
      geofac(0)=geofac(1)
      georotrate(0)=georotrate(1)
      geoalpha(0)=geoalpha(1)
      shat_exp(jmaxm)=shat_exp(jmaxm-1)
      georotrate(jmaxm)=georotrate(jmaxm-1)
      geoalpha(jmaxm)=geoalpha(jmaxm-1)
      geofac(jmaxm)=geofac(jmaxm-1)
      drhodr(jmaxm)=drhodr(jmaxm-1)
      drhodrrrho(jmaxm)=drhodrrrho(jmaxm-1)
c
c      do j=0,jmaxm
c        write(*,60) j,rho(j),rmin_exp(j),vol_exp(j),
c     >              gradrho_exp(j),gradrhosq_exp(j), bteff_exp(j)
c      enddo
c      do j=0,jmaxm
c        write(*,60) j,rho(j),drhodr(j),elong_exp(j),
c     >              delta_exp(j),sfarea_exp(j), vol_exp(j)
c      enddo
c      do j=0,jmaxm
c        write(*,60) j,rho(j),ne_exp(j), ni_exp(j), nz_exp(j),
c     >              zeff_exp(j)
c      enddo
c
c... Write out some 0d quantities
c
      if (i_proc.eq.0) then
        write(6,'(a32,2F10.3)') ' Total integrated Pnbe,Pnbi [MW]:',
     >   pbescale*powe_beam_exp(nj_d-1),pbiscale*powi_beam_exp(nj_d-1)
        write(6,'(a32,2F10.3)') ' Total integrated Prfe,Prfi [MW]:',
     >   prfscale*powe_rf_exp(nj_d-1),prfscale*powi_rf_exp(nj_d-1)
        write(6,'(a32,2F10.3)') ' Total integrated Pcx,Pion [MW]: ',
     >   powi_cx_exp(nj_d-1),powi_ion_exp(nj_d-1)
        write(6,'(a32,2F10.3)') ' Total integrated Poh,Prad [MW]: ',
     >   powe_oh_exp(nj_d-1),xrad_exp*powe_rad_exp(nj_d-1)
        write(6,'(a32,2F10.3)') ' Total integrated Pexch-e,i [MW]:',
     >   -pow_ei_exp(nj_d-1),pow_ei_exp(nj_d-1)
        write(6,'(a32,2F10.3)') ' Total integrated Pwdot-e,i [MW]:',
     >   (1.D0-xwdot)*powe_wdot_exp(nj_d-1),
     >   (1.D0-xwdot)*powi_wdot_exp(nj_d-1)
        if(ialpha.eq.1) then
          write(6,'(a31,1x,2F10.3)') 'Total integrated Pfus-e,i [MW]:',
     >     powe_fus_exp(nj_d-1),powi_fus_exp(nj_d-1)
        endif
        write(6,'(a32,2F10.3)') ' Total integrated Powe,Powi [MW]: ',
     >   powe_exp(nj_d-1),powi_exp(nj_d-1)
        write(6,'(a32,2F10.3)') ' Total integrated neutral source: ',
     >   flow_exp(nj_d-1)
        write(6,'(a32,2F10.3)') ' Line-average electron density: ',
     >   nebar_exp
      endif
c
c..Diagnostic printout of power flows to external file
c
      if (lprint_pflow .eq. 1) then
      open (40,file='pflow.out',status='unknown')
      write(40,90) nj_d-1
      do j=1,nj_d-1
        write(40,100) j, rho(j), powi_beam_exp(j), -powi_ion_exp(j),
     &                powi_cx_exp(j), pow_ei_exp(j), 
     &                -powi_wdot_exp(j), powi_exp(j)
      enddo
      write(40,55) nj_d-1
      do j=1,nj_d-1
        write(40,100) j, rho(j), powe_beam_exp(j), powe_oh_exp(j),
     &                -powe_ion_exp(j), -powe_rad_exp(j),  
     &                pow_ei_exp(j), -powe_wdot_exp(j), powe_exp(j)
      enddo
      close(40)
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 50   format(i2,2x,0p1f4.2,0p6e13.5)
 55   format(i2,2x,0p6e13.5)
 60   format(i2,2x,0p1f6.4,1p6e13.5)
 80   format(' Warning: grad-rho missing, computing it ...')
 90   format(i2,2x,'rho',5x,'powibeam',5x,'powiion',6x,'powicx',7x,
     &       'powei',8x,'powiwdot',5x,'powitot')
 100  format(i2,2x,0p1f4.2,1p8e13.5)
c
      end

