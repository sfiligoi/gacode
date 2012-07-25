c@gridsetup.f
c jek 18-Jan-11
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c... setup transport grid
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine gridsetup(mxgrid)
c
      use legacyread_interface
c
      implicit none
c      include 'mpif.h'
c      include '../inc/tport.m'
cc      include 'data_d.m'
cc      include 'data_exp.m'
c      include '../inc/model.m'
cc      include 'input.m'
cc      include 'glf.m'
c      include '../inc/ptor.m'
c
      integer j, aj, nex, ngg, mxgrid
      real*8 rhox(0:nj-1)
cc      real*8 rho(0:nj-1)
      real, allocatable, dimension(:) :: r_new,
     &   fcap_new, gcap_new, hcap_new, rbp_new, bp0_new,
     &   te_new, ti_new, q_new, ene_new,
     &   en1_new, en2_new, sion_new, 
     &   srecom_new, scx_new, sbcx_new,
     &   s_new, dudtsv_new, enbeam_new,
     &   enn1_new, ennw1_new, ennv1_new,
     &   sbion_new, sbeam_new, 
     &   curden_new, zeff_new, 
     &   angrot_new, chieinv_new, chiinv_new, 
     &   dpedtc_new, dpidtc_new, 
     &   qconde_new, qcondi_new,
     &   qconve_new, qconvi_new, qbeame_new,
     &   qbeami_new, qdelt_new, qrad_new,
     &   qohm_new, qrfe_new, qrfi_new, 
     &   qione_new, qioni_new, qcx_new, 
     &   qfuse_new, qfusi_new, 
     &   rmaj_new, rmin_new, 
     &   psivolp_new, elongx_new, deltax_new, 
     &   sfareanpsi_new, grho1_new,
     &   grho2_new, xkineo_new,
     &   xb2_new, xbm2_new, xngrth_new,
     &   xgrbm2_new, fm1_new, fm2_new,
     &   fm3_new, fhat_new, torque_new,
     &   ptot_new, pfast_new, er_new
      allocate (r_new(nj),
     &   fcap_new(nj), gcap_new(nj), hcap_new(nj), rbp_new(nj),
     &   bp0_new(nj), te_new(nj), ti_new(nj), q_new(nj), ene_new(nj),
     &   en1_new(nj), en2_new(nj), sion_new(nj), 
     &   srecom_new(nj), scx_new(nj), sbcx_new(nj),
     &   s_new(nj), dudtsv_new(nj), enbeam_new(nj),
     &   enn1_new(nj), ennw1_new(nj), ennv1_new(nj),
     &   sbion_new(nj), sbeam_new(nj), 
     &   curden_new(nj), zeff_new(nj), 
     &   angrot_new(nj), chieinv_new(nj), chiinv_new(nj), 
     &   dpedtc_new(nj), dpidtc_new(nj), 
     &   qconde_new(nj), qcondi_new(nj),
     &   qconve_new(nj), qconvi_new(nj), qbeame_new(nj),
     &   qbeami_new(nj), qdelt_new(nj), qrad_new(nj),
     &   qohm_new(nj), qrfe_new(nj), qrfi_new(nj), 
     &   qione_new(nj), qioni_new(nj), qcx_new(nj), 
     &   qfuse_new(nj), qfusi_new(nj), 
     &   rmaj_new(nj), rmin_new(nj), 
     &   psivolp_new(nj), elongx_new(nj), deltax_new(nj), 
     &   sfareanpsi_new(nj), grho1_new(nj),
     &   grho2_new(nj), xkineo_new(nj),
     &   xb2_new(nj), xbm2_new(nj), xngrth_new(nj),
     &   xgrbm2_new(nj), fm1_new(nj), fm2_new(nj),
     &   fm3_new(nj), fhat_new(nj), torque_new(nj),
     &   ptot_new(nj), pfast_new(nj), er_new(nj) )
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c      if(i_proc.eq.0) then
         write(*,10) ' nj_d = ',nj_d
        if(nj_d.gt.nj) write(*,*) 'Warning: nj_d > nj'
c      endif
cc      jmaxm=mxgrid
c
c... experimental rho_hat grid
c
      do j=0,nj
        rhox(j)=0.D0
      enddo
c
      do j=1,nj_d
        rhox(j)=REAL(j-1)/REAL(nj_d-1)
      enddo
cc      rhox(0)=1.D-6
cc      rhox(nj_d)=rhox(nj_d)-1.D-6
c
c... setup target grid
      do j=1,mxgrid+1
        rho_d(j)=REAL(j-1)/REAL(mxgrid)
      enddo
cc      rho_d(1)=1.D-6
cc      rho_d(mxgrid+1)=rho_d(mxgrid+1)-1.D-6
c
c... Interpolate profiles to new grid if 
c    mxgrid > nj_d-1
c
      if (mxgrid .gt.(nj_d-1)) then
cc        write(*,25)
        nex=mxgrid+1
c        do j=0,51
c          write(*,'(i3,3f8.4)') j, rho_d(j), rhox(j), r_d(j)
c        enddo
c
        ngg=nj_d
c       call w_lin_interp_r8(ngg,rhox,r_d,nex,rho_d,r_new,iflag,msg)
        call inter_cspl(ngg,rhox,r_d,nex,rho_d,r_new)
        call inter_cspl(ngg,rhox,fcap_d,nex,rho_d,fcap_new)
        call inter_cspl(ngg,rhox,gcap_d,nex,rho_d,gcap_new)
        call inter_cspl(ngg,rhox,hcap_d,nex,rho_d,hcap_new)
        call inter_cspl(ngg,rhox,rbp_d,nex,rho_d,rbp_new)
        call inter_cspl(ngg,rhox,bp0_d,nex,rho_d,bp0_new)
        call inter_cspl(ngg,rhox,te_d,nex,rho_d,te_new)
        call inter_cspl(ngg,rhox,ti_d,nex,rho_d,ti_new)
        call inter_cspl(ngg,rhox,q_d,nex,rho_d,q_new)
        call inter_cspl(ngg,rhox,ene_d,nex,rho_d,ene_new)
        call inter_cspl(ngg,rhox,en_d(1:j,1),nex,rho_d,en1_new)
        call inter_cspl(ngg,rhox,en_d(1:j,2),nex,rho_d,en2_new)
        call inter_cspl(ngg,rhox,sion_d(1:j,1),nex,rho_d,sion_new)
        call inter_cspl(ngg,rhox,srecom_d(1:j,1),nex,rho_d,srecom_new)
        call inter_cspl(ngg,rhox,sbcx_d(1:j,1),nex,rho_d,sbcx_new)
        call inter_cspl(ngg,rhox,scx_d(1:j,1),nex,rho_d,scx_new)
        call inter_cspl(ngg,rhox,s_d(1:j,1),nex,rho_d,s_new)
        call inter_cspl(ngg,rhox,dudtsv_d(1:j,1),nex,rho_d,dudtsv_new)
        call inter_cspl(ngg,rhox,enbeam_d,nex,rho_d,enbeam_new)
        call inter_cspl(ngg,rhox,enn_d(1:j,1),nex,rho_d,enn1_new)
        call inter_cspl(ngg,rhox,ennw_d(1:j,1),nex,rho_d,ennw1_new)
        call inter_cspl(ngg,rhox,ennv_d(1:j,1),nex,rho_d,ennv1_new)
        call inter_cspl(ngg,rhox,sbion_d,nex,rho_d,sbion_new)
        call inter_cspl(ngg,rhox,sbeam_d,nex,rho_d,sbeam_new)
        call inter_cspl(ngg,rhox,curden_d,nex,rho_d,curden_new)
        call inter_cspl(ngg,rhox,zeff_d,nex,rho_d,zeff_new)
        call inter_cspl(ngg,rhox,angrot_d,nex,rho_d,angrot_new)
        call inter_cspl(ngg,rhox,chieinv_d,nex,rho_d,chieinv_new)
        call inter_cspl(ngg,rhox,chiinv_d,nex,rho_d,chiinv_new)
        call inter_cspl(ngg,rhox,dpedtc_d,nex,rho_d,dpedtc_new)
        call inter_cspl(ngg,rhox,dpidtc_d,nex,rho_d,dpidtc_new)
        call inter_cspl(ngg,rhox,qconde_d,nex,rho_d,qconde_new)
        call inter_cspl(ngg,rhox,qcondi_d,nex,rho_d,qcondi_new)
        call inter_cspl(ngg,rhox,qconve_d,nex,rho_d,qconve_new)
        call inter_cspl(ngg,rhox,qconvi_d,nex,rho_d,qconvi_new)
        call inter_cspl(ngg,rhox,qbeame_d,nex,rho_d,qbeame_new)
        call inter_cspl(ngg,rhox,qbeami_d,nex,rho_d,qbeami_new)
        call inter_cspl(ngg,rhox,qdelt_d,nex,rho_d,qdelt_new)
        call inter_cspl(ngg,rhox,qrfe_d,nex,rho_d,qrfe_new)
        call inter_cspl(ngg,rhox,qrfi_d,nex,rho_d,qrfi_new)
        call inter_cspl(ngg,rhox,qione_d,nex,rho_d,qione_new)
        call inter_cspl(ngg,rhox,qioni_d,nex,rho_d,qioni_new)
        call inter_cspl(ngg,rhox,qcx_d,nex,rho_d,qcx_new)
        call inter_cspl(ngg,rhox,qfuse_d,nex,rho_d,qfuse_new)
        call inter_cspl(ngg,rhox,qfusi_d,nex,rho_d,qfusi_new)
        call inter_cspl(ngg,rhox,qrad_d,nex,rho_d,qrad_new)
        call inter_cspl(ngg,rhox,qohm_d,nex,rho_d,qohm_new)
        call inter_cspl(ngg,rhox,rmajavnpsi_d,nex,rho_d,rmaj_new)
        call inter_cspl(ngg,rhox,rminavnpsi_d,nex,rho_d,rmin_new)
        call inter_cspl(ngg,rhox,psivolp_d,nex,rho_d,psivolp_new)
        call inter_cspl(ngg,rhox,elongx_d,nex,rho_d,elongx_new)
        call inter_cspl(ngg,rhox,deltax_d,nex,rho_d,deltax_new)
        call inter_cspl(ngg,rhox,sfareanpsi_d,nex,rho_d,sfareanpsi_new)
        call inter_cspl(ngg,rhox,grho1npsi_d,nex,rho_d,grho1_new)
        call inter_cspl(ngg,rhox,grho2npsi_d,nex,rho_d,grho2_new)
        call inter_cspl(ngg,rhox,torque_d,nex,rho_d,torque_new)
        call inter_cspl(ngg,rhox,ptot_d,nex,rho_d,ptot_new)
        call inter_cspl(ngg,rhox,pfast_d,nex,rho_d,pfast_new)
cc        if(ncl_flag.ne.0)then
          call inter_cspl(ngg,rhox,xb2_d,nex,rho_d,xb2_new)
          call inter_cspl(ngg,rhox,xbm2_d,nex,rho_d,xbm2_new)
          call inter_cspl(ngg,rhox,xngrth_d,nex,rho_d,xngrth_new)
          call inter_cspl(ngg,rhox,xgrbm2_d,nex,rho_d,xgrbm2_new)
          call inter_cspl(ngg,rhox,fm1_d,nex,rho_d,fm1_new)
          call inter_cspl(ngg,rhox,fm2_d,nex,rho_d,fm2_new)
          call inter_cspl(ngg,rhox,fm3_d,nex,rho_d,fm3_new)
          call inter_cspl(ngg,rhox,fhat_d,nex,rho_d,fhat_new)
cc        endif
c
c... map interpolated variables back to data (_d) variables
c
c        do j=1,jmaxm+1
c          write(*,50) j, rhox(j-1), r_d(j), ti_d(j), 
c     &                rho_d(j), r_new(j), ti_new(j)
c        enddo
        do j=1,nex
          r_d(j)=r_new(j)
          fcap_d(j)=fcap_new(j)
          gcap_d(j)=gcap_new(j)
          hcap_d(j)=hcap_new(j)
          rbp_d(j)=rbp_new(j)
          bp0_d(j)=bp0_new(j)
          te_d(j)=te_new(j)
          ti_d(j)=ti_new(j)
          q_d(j)=q_new(j)
          ene_d(j)=ene_new(j)
          en_d(j,1)=en1_new(j)
          en_d(j,2)=en2_new(j)
          sion_d(j,1)=sion_new(j)
          srecom_d(j,1)=srecom_new(j)
          sbcx_d(j,1)=sbcx_new(j)
          scx_d(j,1)=scx_new(j)
          s_d(j,1)=s_new(j)
          dudtsv_d(j,1)=dudtsv_new(j)
          enbeam_d(j)=enbeam_new(j)
          enn_d(j,1)=enn1_new(j)
          ennw_d(j,1)=ennw1_new(j)
          ennv_d(j,1)=ennv1_new(j)
          sbion_d(j)=sbion_new(j)
          sbeam_d(j)=sbeam_new(j)
          curden_d(j)=curden_new(j)
          zeff_d(j)=zeff_new(j)
          angrot_d(j)=angrot_new(j)
          chieinv_d(j)=chieinv_new(j)
          chiinv_d(j)=chiinv_new(j)
          dpedtc_d(j)=dpedtc_new(j)
          dpidtc_d(j)=dpidtc_new(j)
          qconde_d(j)=qconde_new(j)
          qcondi_d(j)=qcondi_new(j)
          qconve_d(j)=qconve_new(j)
          qconvi_d(j)=qconvi_new(j)
          qbeame_d(j)=qbeame_new(j)
          qbeami_d(j)=qbeami_new(j)
          qdelt_d(j)=qdelt_new(j)
          qrfe_d(j)=qrfe_new(j)
          qrfi_d(j)=qrfi_new(j)
          qione_d(j)=qione_new(j)
          qioni_d(j)=qioni_new(j)
          qcx_d(j)=qcx_new(j)
          qfuse_d(j)=qfuse_new(j)
          qfusi_d(j)=qfusi_new(j)
          qrad_d(j)=qrad_new(j)
          qohm_d(j)=qohm_new(j)
          rmajavnpsi_d(j)=rmaj_new(j)
          rminavnpsi_d(j)=rmin_new(j)
          psivolp_d(j)=psivolp_new(j)
          elongx_d(j)=elongx_new(j)
          deltax_d(j)=deltax_new(j)
          sfareanpsi_d(j)=sfareanpsi_new(j)
          grho1npsi_d(j)=grho1_new(j)
          grho2npsi_d(j)=grho2_new(j)
          torque_d(j)=torque_new(j)
          ptot_d(j)=ptot_new(j)
          pfast_d(j)=pfast_new(j)
c          er_d(j)=er_new(j)
        enddo
cc        if (ncl_flag .eq. 1) then
          do j=1,nex
            xb2_d(j)=xb2_new(j)
            xbm2_d(j)=xbm2_new(j)
            xngrth_d(j)=xngrth_new(j)
            xgrbm2_d(j)=xgrbm2_new(j)
            fm1_d(j)=fm1_new(j)
            fm2_d(j)=fm2_new(j)
            fm3_d(j)=fm3_new(j)
            fhat_d(j)=fhat_new(j)
          enddo
cc        endif
c
c deallocate temp arrays
c
      deallocate( r_new, 
     &   fcap_new, gcap_new, hcap_new, rbp_new, bp0_new,
     &   te_new, ti_new, q_new, ene_new,
     &   en1_new, en2_new, sion_new, 
     &   srecom_new, scx_new, sbcx_new,
     &   s_new, dudtsv_new, enbeam_new,
     &   enn1_new, ennw1_new, ennv1_new,
     &   sbion_new, sbeam_new, 
     &   curden_new, zeff_new, 
     &   angrot_new, chieinv_new, chiinv_new, 
     &   dpedtc_new, dpidtc_new, 
     &   qconde_new, qcondi_new,
     &   qconve_new, qconvi_new, qbeame_new,
     &   qbeami_new, qdelt_new, qrad_new,
     &   qohm_new, qrfe_new, qrfi_new, 
     &   qione_new, qioni_new, qcx_new, 
     &   qfuse_new, qfusi_new, 
     &   rmaj_new, rmin_new, 
     &   psivolp_new, elongx_new, deltax_new, 
     &   sfareanpsi_new, grho1_new,
     &   grho2_new, xkineo_new,
     &   xb2_new, xbm2_new, xngrth_new,
     &   xgrbm2_new, fm1_new, fm2_new,
     &   fm3_new, fhat_new, torque_new,
     &   ptot_new, pfast_new, er_new)
c
cjek nj_d -> mxgrid+1 for arbitrary grid size
c
      nj_d=mxgrid+1
c
      endif
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 10   format(a9,i3)             ! common write/read  format
 25   format(' Interpolating ONETWO iterdb grid to new grid ...')
 50   format(i2,2x,0p1f4.2,0p6f10.5)
c
      end
