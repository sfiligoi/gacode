      subroutine init
c
c initialize variables
c
      implicit none
c
      include '../inc/input.m'
      include '../inc/tport.m'
      include '../inc/model.m'
      include '../inc/data.m'
      include '../inc/ptor.m'
      include '../inc/vtrans.m'
      include '../inc/glf.m'
c
      integer j, k
c
      do j=0,jmaxmt
        powe_exp(j)=0.D0
        powi_exp(j)=0.D0
        pow_ei_exp(j)=0.D0
      enddo
      do j=1,mxgrd
        Peaux(j)=0.D0
        Piaux(j)=0.D0
      enddo
      do j=1,mxfds
        do k=1,mxgd
          Told(j,k)=0.D0
          Tnew(j,k)=0.D0
          vrho(j,k)=0.D0
          S(j,k)=0.D0
          grow(j,k)=0.D0
          CONV(j,k)=0.D0
        enddo
      enddo
      do j=1,mxfields
        sramp_pt(j)=0.D0
        s0_pt(j)=1.D0
        s1_pt(j)=1.D0
        smult(j)=1.D0
        ave_field(j)=0.D0
      enddo
      do j=1,8
        upsrc(j)=0
      enddo
      amassgas_exp=2.0
      amassimp_exp=12.0
      idatzero=0
      iv_dk=2
      alpha_e_dk=0.D0
      alpha_star_dk=0.D0
      time_upd=99.D0
      endtime_pt=0.D0
      restart_pt=-1.D0
      time_series = 0
      ipramp = 0
      ismooth = 0
      ismooth_all=0
      ilog=0
      ifixeta=0
      dvflag=1
      ave_dv = 0.D0
      adiff_dv = 0.D0
      d_art = 0.D0
      dt_implicit=1.D0
      zeff_e = 0.D0
      ne_gain_pt=0.0
      te_gain_pt=0.0
      ti_gain_pt=0.0
      wall_mult=1.0
      te_bc_mult=1.0
      ti_bc_mult=1.0
c      if(iexch.eq.0) iexch=1
      irot1=1
      irot2=1
      iohm=0
      irad=0
      xconv=0.D0
      xwdot=0.D0
      xsdot=0.D0
      sigma=0.D0
      ichiep=0
      pimpa=12.D0
      pimpz=6.D0
      echconv=1.D0
      ipptot=0 ! don't include fast ions in ptot in rep_iter.f
      itorque=0
      iexb=0
      igyro=0
      istringer_nc=0
      ismoo_ne=0.D0
      ismoo_ni=0.D0
      ismoo_nf=0.D0
      ismoo_zeff=0.D0
      ismoo_q=0.D0
      ismoo_vrot=0.D0
      smoo_tim=0.D0
      smoo_vphim=0.D0
      istep_smoo=9999
      ifilter = 0
      xion_exp=1.D0
      xrad_exp=1.D0
      didledge_n=0.D0
      didledge_te=0.D0
      didledge_ti=0.D0
      didledge_vphi=0.D0
      vexb_bc=0.D0
      time1_pow=0.D0
      time2_pow=99.D0
      time3_pow=99.D0
      pbescale2=1.D0
      pbiscale2=1.D0
      pbescale3=1.D0
      pbiscale3=1.D0
      dilution_model=0
c
      return
      end
