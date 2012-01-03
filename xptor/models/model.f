c@model.f J. Kinsey
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  Transport models                                           8/14/01
c
c  This routine groups together the various turbulent transport
c  models. The model is chosen using the variable 'imodel'.
c
c      imodel          description
c        0          simple critical gradient model
c        2          IFS/PPPL ITG model
c        3          Itoh-Itoh-Fukuyama current diffusive model
c        4          Culham semi-empirical model
c        6          Multi-mode 1995 model
c        7          Mixed-shear semi-empirical model
c        8          old GLF23 v1.0 ITG/TEM/ETG model
c        81         updated GLF23 v1.0 ITG/TEM/ETG model
c        99         chii=chie=1.0 m^2/s
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c Changes :
c
c 8/4/00 jek added grid scale diffusion
c 5/8/00 jek added Culham model
c 4/24/00 jek added IIF model
c 1/17/00 jek pass corot*vphi_m to callglf2d.f instead of just vphi_m
cmnt   given gradients model returns flows at a given radius       
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine model
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
      integer nerr, i, ik, j, k, jshoot, jmm, igrad, idvflag
      integer k_order, k_potato, m_i, m_z, im, iza, iz, idum(8)
      real*8 lnlam, gfac
      real*8 ve(0:jmaxm),vpar(0:jmaxm),vper(0:jmaxm)
      real*8 chiemi_ms, ceqei_ms, facbr_ms, facgb_ms
      real*8 fshea_ms, dbohm_ms, gyro_ms, rlte_ms, rlti_ms
      real*8 chienem, chietim, chiitem, chiinem, difftem, difftim
      real*8 zepsne, zepsni, zepste, zepsti, zgne, zgni, zgte, zgti
      real*8 zpmne_wn, epsnein_wn, gnein_wn, gnhin_wn, gnzin_wn
      real*8 gniin_wn, gnsin_wn, gnqin_wn
      real*8 gtein_wn, gthin_wn, gtzin_wn
      real*8 zgyrfi_wn, cs_wn, nue, gamma_test_max, zcc_kb
      real*8 zcmu0_kb, zceps0_kb, zepslon_kb, sqrtelong_kb
      real*8 zprth, zdprth, zsgdpr, zdpdr, cthery3_kb, zlgeps_kb
      real*8 betaein_kb, betaiin_kb, zbeta_kb, qin_kb, shearin_k
      real*8 zbprim_kb, zbcoef1_kb, zbcoef2_kb, zbc1_kb, zbc2_kb
      real*8 zgmax_kb, zgpr_kb
      real*8 zbpbc1, zfbthn_kb, zdk_kb, zpmne_q, zpmni_q
      real*8 aiwt_jp1, xnimp_jp1, xnimp, vstar_sign, diffgb_local
      real*8 exchgb_local, alpha_e, x_alpha, exchm, fac_imp_flow
      real*8 zpmnix, thresitg
      real*8 c_den, c_potb, c_potl, p_b2, p_bm2, p_eb, p_fhat
      real*8 p_fm(3), p_ft, p_grbm2, p_ngrth, p_grphi, p_gr2phi
      real*8 elemas, promas, drho, art_diff
c      real*8 radfunc(0:jmaxmt)
c
c IFS variables
      real*8 RLT_dk, RLN_dk, RLNe_dk, q_dk, akappa_dk, shat_dk,
     &       zth_dk, xnbeam_dk, tau_dk, eps_dk, xnu_dk, rho_i_dk,
     &       v_ti_dk, rmajor_dk, g_perp_dk, RLTcrit_dk, RLTcrit2_dk,
     &       chi_0_dk, g_dk, chi_i_dk, chi_e_dk, RLTrotshear_dk,
     &       gamma_dk, rotstab,rotshear_dk
      real*8 vti_dk, omegaci_dk, rhoi_dk, xnu_ntcc_dk
c
c IIF variables
      real*8 v_alfven_ii, c_skindepth_ii, alpha_ii, akap_ii,
     &       gamma_e_ii, h_ii, shat_p, shat_ii, f_ii, g1_ii, chi_ii
c
      real*8 z_j7kv
c
c Culham variables
      real*8 tau_cul, amin_cul, eps_cul, ni_cul, ni1_cul,
     &       p_cul, p1_cul, grdte_cul, grdti_cul, grdp_cul,
     &       rlti_cul, rlte_cul, rlp_cul, alpha_cul, falpha_cul,
     &       chie1_cul, chie2_cul, chii_cul
c
c MMM95 variables
      integer npoints, nnout, lsuper, lreset, lprint_wn
      integer jr, lswitch(8)
      real*8 cswitch(25)
      real*8 nitotm, fig_mm(4), fkb_mm(4), frb_mm(4)
      real*8 rmin_mm, rmaj_mm, elong_mm, ne_mm, nh_mm, nz_mm, ns_mm,
     &       nse_mm, nitot_mm, zeff_mm, te_mm, ti_mm, q_mm, btor_mm,
     &       zimpz_mm, zimpa_mm, zgasa_mm, zgasz_mm, zaimass_mm
c
      real*8 zthiig,   zthdig,   ztheig,   zthzig
     & , zthirb,   zthdrb,   ztherb,   zthzrb
     & , zthikb,   zthdkb,   zthekb,   zthzkb
      real*8 zgamma_wn(12), zomega_wn(12)
      real*8 shear_wn, zscyl_wn
      real*8 ztau_nco, zrho_nco
c
c
c FORCEBAL variables
      character*15 device
      integer mxnr_r, nr_r
      parameter(mxnr_r=130)
      include '../inc/pamx_ms.inc'
      include '../inc/pamx_mz.inc'
      include '../inc/comfbl.inc'
      real rhot_r(mxnr_r), time_ncl
      integer mx_ni, xneo
      parameter(mx_ni=5)
      character*3 ci0(mx_ni), cim1, cim2
c
      save idvflag
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
c   kev/sec per MW       
c       kevdsecpmw=1.6022D-19*1.0D3*1.D-6
c
c   electron and proton masses
c
        elemas=9.1095e-31
        promas = 1.6726e-27
        z_j7kv=1.6022e-16
c
c   unit of diffusion in m**2/sec   (computed in derived modprofiles)
c       cgyrobohm_m(jm)=1.D-4*
c    >  979.D3*(tem*1.D3)**.5D0/(arho_exp*100.)
c    >  *(102.D0*(tem*1.D3)**.5D0/bt_exp/1.D4)**2*(amassgas_exp)**.5D0
c   unit of density is 10**13 cm**-3 or 10**19 m**-3
c   unit of arho_exp and rmajor_exp is meters
c   sfactor(j)=
c     >      2.*pi_m*arho_exp*rho(j)*h_exp(j)*2.*pi_m*rmajor_exp 
c   units of meter**2
c   unit of power is MW
c   unit of flow  is MW/kev=kamp
c
c...initialize diffusivities
c
        diffnem=0.D0
        chietem=0.D0
        chiitim=0.D0
        etaphim=0.D0
c
 9000 format (1i6,7e10.3)
c       
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cmnt                         the model
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cmnt
cmnt supply model for chietem,chietim,chienem
cmnt                  chiitem,chiitim,chiinem
cmnt                  difftem,difftim,diffnem
cmnt     
cmnt    chi_s and diff_s must be in meters**2/sec units
cmnt     and recall chi_s refer to total energy flow
cmnt
cmnt    if the model chi_s refer to "heat conduction" flow
cmnt    then a convection term xconv*1.5D0*t_m*flow_exp is added.
cmnt    normally input xconv=0. otherwise xconv=1. or 5./3.
cmnt    
cmnt    it is also possible to build convection into the model
cmnt    with "aconv".  aconv and xconv should not be double counted.
cmnt 
cmnt note: can use diagonal forms with off diagonal dependence on
cmnt zpmte,zpmti,zpmne intrinsic to the diagonal elements as in sample
cmnt normall models are written in diagonal for with dependence on off 
cmnt diagonal gradients implicit
cmnt
cmnt note: when flow is large anomalous D-i exchange should be added
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (imodel.eq.0) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c simple model
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      thresitg= (zpmti-zpticrit*((1.D0-amod0(1))+amod0(1)*tim/tem))
      if (thresitg.gt.0.)  thresitg=thresitg**epsthres
      if (thresitg.le.0.)  thresitg=thresitg0
      if (epsthres.lt.0.)  thresitg=1.D0
c
c radial profile function
c
      betarad=1.D-3
      alfarad=4.D0
      epsrad=20.D0
      do j=0,jmaxm
        radfunc(j)=(1.D0+epsrad*rho(j)**alfarad*
     &  (1.D0-rho(j)**2.D0)**betarad)/(1.D0+epsrad)
      enddo
c
      diffnem=cmodel*(dntem+dttem*zpmte/zpmne+
     > (ctitge+cnitge*zpmne/zpmte)*thresitg)*
     >   radfunc(jm)*cgyrobohm_m(jm)
      difftem=0.D0
      difftim=0.D0
c     
      chietem=cmodel*(ctteme+cnteme*zpmne/zpmte 
     > +(ctitge+cnitge*zpmne/zpmte)*thresitg)*
     >   radfunc(jm)*cgyrobohm_m(jm)
     > +3.D0/2.D0*aconv*diffnem*zpmne/zpmte
      chietim=0.D0
      chienem=0.D0
c
      chiitem=0.D0
      chiitim=cmodel*(cttemi+cntemi*zpmne/zpmti
     > +(ctitgi+cnitgi*zpmne/zpmti)*thresitg)*
     >  radfunc(jm)*cgyrobohm_m(jm)
     > +3.D0/2.D0*aconv*diffnem*zpmne/zpmti
      chiinem=0.D0

      etaphim=etaphi0*chiitim
      etaparm=etapar0*chiitim
      etaperm=etaper0*chiitim
c
      if(cmodel.lt.0) then
       chietem=-cmodel
       chiitim=-cmodel
       if(cneo.ne.0.and.ineo.eq.-20) chiitim=0.D0
      endif
c
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
       etagb_phi_m(jm)=etaphim/cgyrobohm_m(jm)
       etagb_par_m(jm)=etaparm/cgyrobohm_m(jm)
       etagb_per_m(jm)=etaperm/cgyrobohm_m(jm)
c
       chie_e_gb_m(jm)=chie_e_gf*cmodel*gfac
c
       diff_m(jm)=diffnem
       chie_m(jm)=chietem
       chii_m(jm)=chiitim
       etapar_m(jm)=etaparm
       etaper_m(jm)=etaperm
       etaphi_m(jm)=etaphim  
c
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (imodel.eq.2) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    IFS PPPL (Dorland Kotchenreuther) model
c    Two versions are implemented - The 1994 original version (ip_chi1)
c    and the 1996 Varenna version (ip_chi2) with both velocity and 
c    diamagnetic shear stabilization plus an additional elongation factor.
c
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
c.. use NTCC switches if itest_ntcc > 0
c   Puts limits on rlne, q, shat, nu, tau, eps
c
       if(itest_ntcc.gt.0) then
         if(q_exp(jm).lt.0.7D0) q_dk=0.7D0
         if(q_exp(jm).gt.8.0D0) q_dk=8.0D0
         if(shat_dk.lt.0.5D0) shat_dk=0.5D0
         if(shat_dk.gt.2.0D0) shat_dk=2.0D0
         if(eps_dk.lt.0.1D0) eps_dk=0.1D0
         if(eps_dk.gt.0.3D0) eps_dk=0.3D0
         if(RLN_dk.lt.1.D-6) RLNe_dk=1.D-6
         if(RLN_dk.gt.6.D0) RLNe_dk=6.D0
         vti_dk=dsqrt(tim*1.6022D-16/3.3436D-27)
         omegaci_dk=1.6022D-19*bt_exp/3.3436D-27
         rhoi_dk=vti_dk/omegaci_dk
         if(tim.lt.1.D-4 .or. tem.lt.1.D-4) then
           xnu_ntcc_dk=10.D0
         else
           xnu_ntcc_dk=2.5D0*rmajor_exp*nem/(tem**1.5D0*tim**0.5D0)
           if(xnu_ntcc_dk.lt.0.5D0) xnu_ntcc_dk=0.5D0
           if(xnu_ntcc_dk.gt.10.D0) xnu_ntcc_dk=10.D0
         endif
c        write(*,50) jm, rho(jm), v_ti_dk, vti_dk, rho_i_dk, rhoi_dk
c        write(*,50) jm, rho(jm), xnu_dk, xnu_ntcc_dk
         rho_i_dk=rhoi_dk
         xnu_dk=xnu_ntcc_dk
         if(tem.lt.1.D-4) then
           tau_dk=4.D0
         else
           tau_dk=tim/tem
           if(tau_dk.lt.0.5D0) tau_dk=0.5D0
           if(tau_dk.gt.4.0D0) tau_dk=4.0D0
         endif
       endif
c
c
c..jek 1994 version of IFS model
c
       if(iv_dk.eq.1) then
       call ip_chi1(itest_dk,RLT_dk,RLN_dk,q_dk,
     .  shat_dk,zth_dk,
     .  xnbeam_dk,tau_dk,eps_dk,xnu_dk,
     .  rmajor_dk,rho_i_dk,v_ti_dk,RLTcrit_dk,
     .  RLTcrit2_dk,chi_0_dk,g_dk,chi_i_dk,chi_e_dk)
c
       if (alpha_e_dk.ne.0 .or. alpha_star_dk.ne.0) then
         RLTrotshear_dk=
     >    (alpha_e_dk*abs(egamma_exp(jm))
     >    +alpha_star_dk*abs(vstarp_exp(jm)))
     >    /(0.1D0*arho_exp/rmajor_exp)
c
         if(irotstab.gt.0) then
         RLTrotshear_dk=
     >    (alpha_e_dk*abs(egamma_m(jm))
     >    +alpha_star_dk*abs(vstarp_m(jm)))
     >    /(0.1D0*arho_exp/rmajor_exp)
         endif 
       else
         RLTrotshear_dk=0.D0
       endif
       gamma_dk=0.1D0*arho_exp/rmaj_exp(jm)*
     >       (RLT_dk-RLTcrit_dk+RLTrotshear_dk)
       rotstab=1.D0
       endif
c
c..jek Varenna96 version of IFS model 
c      set shearing rate g_perp zero and compute
c      rotstab outside of ip_chi routine

       if(iv_dk.eq.2) then
       call ip_chi2(itest_ntcc,itest_dk,
     .  RLT_dk,RLN_dk,RLNe_dk,q_dk,akappa_dk,
     .  shat_dk,zth_dk,
     .  xnbeam_dk,tau_dk,eps_dk,xnu_dk,g_perp_dk,
     .  rmajor_dk,rho_i_dk,v_ti_dk,RLTcrit_dk,
     .  RLTcrit2_dk,chi_0_dk,g_dk,
     .  gamma_dk,chi_i_dk,chi_e_dk)
c
       gamma_dk=gamma_dk/csda_m(jm)
c
c..rotational stabilization using either the
c  experimental (irotstab=0) or model (irotstab=1)
c  temperature
c
       rotstab=1.D0
       if ( irotstab.eq.0 ) then
         rotshear_dk=alpha_e_dk*egamma_exp(jm)+
     >              alpha_star_dk*vstarp_exp(jm)
       rotstab=(1.D0-abs(rotshear_dk/
     >         (gamma_dk*csda_m(jm)/csda_exp(jm))))
       elseif (irotstab.gt.0) then
         rotshear_dk=alpha_e_dk*egamma_m(jm)+
     >              alpha_star_dk*vstarp_m(jm)
       rotstab=(1.D0-abs(rotshear_dk/gamma_dk))
       else
         rotshear_dk=0.D0
       endif
c
       rotstab_m(jm)=rotstab
       if (rotstab.lt.0.D0) rotstab=1.D-6
       endif
c
c gamma_dk is a description of max linear growth rate from gks code
       if (gamma_dk.le.0.) gamma_dk=1.D-6
c
       chietem=cmodel_e*chi_e_dk*geofac(jm)
       chietim=0.D0
       chienem=0.D0
c
       chiitem=0.D0
       chiitim=cmodel_i*chi_i_dk*geofac(jm)
       chiinem=0.D0
c
       difftem=0.D0
       difftim=0.D0
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
       chietem=chietem*rotstab
       chiitim=chiitim*rotstab
       diffnem=diffnem*rotstab
c
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
       anrate_m(jm)=gamma_dk
       dnrate_m(jm)=0.D0
       dtnrate_m(jm)=0.D0
       anfreq_m(jm)=0.D0
       dnfreq_m(jm)=0.D0
c
       if(istep_glf.eq.lprint_glfstep) then
         if(idvflag.eq.0) then
c           write(*,51) jm, rho(jm), RLT_dk, RLTcrit_dk
         endif
         idvflag=idvflag+1
         if(idvflag.eq.2) idvflag=0
       endif
      endif

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if (imodel.eq.3) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmnt Itoh-Itoh model
cmnt PRL 72 (1994) 1200
c
      v_alfven_ii=1.D-2*2.18D11*(bt_exp*1.D4)/(nim*1.D13)**0.5D0/
     >   (amassgas_exp)**0.5D0   !m/sec
      c_skindepth_ii=1.D-2*5.31D5/(nem*1.D13)**0.5D0   !meters
c fix later ...much later
      zpmni=zpmne
      alpha_ii=betae_m(jm)*q_exp(jm)**2*rmaj_exp(jm)/arho_exp
     >     *drhodr(jm)
     >  *(zpmte+zpmne+(tim/tem)*(zpmti+zpmni))
      shat_ii=shat_exp(jm)*drhodrrrho(jm)
      akap_ii=-rmin_exp(jm)/rmaj_exp(jm)*(1.D0-1.D0/q_exp(jm))
      gamma_e_ii=egamma_exp(jm)
      h_ii=q_exp(jm)*rmaj_exp(jm)/v_alfven_ii*csda_exp(jm)*gamma_e_ii
      
      shat_p=shat_ii-alpha_ii
      if(shat_p.le.0.) then
       f_ii=(1.D0+akap_ii)**(5.D0/2.D0)/
     >      dsqrt(2.D0*(1.D0-2.D0*shat_p)*
     >      (1.D0-2.D0*shat_p+3.D0*shat_p**2*(1.D0+akap_ii)))
      else
       f_ii=(1.D0+akap_ii)**(5.D0/2.D0)*(1.D0+9.D0*dsqrt(2.D0)*
     >      shat_p**(5.D0/2.D0))/dsqrt(2.D0)/(1.D0-2.D0*shat_p+
     >      3.D0*shat_p**2*(1.D0+akap_ii)+2.D0*shat_p**3)
      endif
  
c check this: the gamma_e effect is left out with g1_ii=0
c PRL give coomplicated formula for g1_ii but Fuhayama leaves
c it out 
c   
       g1_ii=0.D0
       f_ii=1.D0/(1.D0+g1_ii*h_ii**2)*f_ii
c 
       chi_ii=f_ii*(abs(alpha_ii))**(3.D0/2.D0)*
     >      (c_skindepth_ii)**2*v_alfven_ii/q_exp(jm)/rmaj_exp(jm)
c      
       chietem=12.D0*cmodel_e*geofac(jm)*cmodel*chi_ii
       chietim=0.D0
       chienem=0.D0
       chiitem=0.D0
       chiitim=12.D0*cmodel_i*geofac(jm)*cmodel*chi_ii
       chiinem=0.D0
       difftem=0.D0
       difftim=0.D0
       diffnem=0.D0 
c  
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)  
c
      endif
c
      if (imodel.eq.4 .or. imodel.eq.41) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjek Culham model
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       tau_cul=ti_m(jm)/max(te_m(jm),1.D-10)
       amin_cul=amin_d*(kappa_d+1.D0)/2.D0
       eps_cul=amin_cul*rho(jm)/rmaj_exp(jm)
       ni_cul=ni_exp(jm)+nz_exp(jm)
       ni1_cul=ni_exp(jm+1)+nz_exp(jm+1)
       p_cul=(ne_m(jm)*te_m(jm)+ni_cul*ti_m(jm))*1.D19*1.D3
       p1_cul=(ne_m(jm+1)*te_m(jm+1)+ni1_cul*ti_m(jm+1))*1.D19*1.D3
       grdte_cul=(te_exp(jm)-te_exp(jm+1))*1.D3 / (rho(jm)-rho(jm+1))
       grdti_cul=(ti_exp(jm)-ti_exp(jm+1))*1.D3 / (rho(jm)-rho(jm+1))
       grdp_cul=(p_cul-p1_cul) / (rho(jm)-rho(jm+1))
       rlti_cul=rmaj_exp(jm)*abs(grdti_cul/(ti_m(jm)*1.D3*amin_cul))
       rlte_cul=rmaj_exp(jm)*abs(grdte_cul/(te_m(jm)*1.D3*amin_cul))
       rlp_cul=rmaj_exp(jm)*abs(grdp_cul/(p_cul*amin_cul))
       alpha_cul=4.03D-25*(p_cul*q_exp(jm)**2*rlp_cul/bt_exp**2)
       if (imodel.eq.4) then
         falpha_cul=alpha_cul**3.1D0/(alpha_cul**2+1.D-2)
       else
         falpha_cul=alpha_cul**3.1D0/(alpha_cul**3+1.D-3)
       endif
c
       chie1_cul=1.186D18*(dsqrt(te_m(jm)*1.D3)*eps_cul*q_exp(jm))/
     &           (ne_m(jm)*1.D19*rmaj_exp(jm)*dsqrt(amassgas_exp))*
     &           (1.D0+rlte_cul/2.D0)
       if (imodel.eq.4) then
         chie2_cul=2.372D-5*((te_m(jm)*1.D3)**(3.D0/2.D0)*q_exp(jm))/
     &           (bt_exp**2*rmaj_exp(jm)*dsqrt(amassgas_exp))*
     &           (eps_cul*q_exp(jm)*rlp_cul*(1.D0+tau_cul)+1.D0)
       else
         chie2_cul=4.75D-5*((te_m(jm)*1.D3)**(3.D0/2.D0)*q_exp(jm))/
     &           (bt_exp**2*rmaj_exp(jm)*dsqrt(amassgas_exp))*
     &           (eps_cul*rlp_cul*(1.D0+tau_cul)+1.D0)
       endif
c
       if (imodel.eq.4) then
         chii_cul=1.D-2*(ti_m(jm)*1.D3*rlti_cul*tau_cul*falpha_cul)/
     &          (dsqrt(amassgas_exp)*bt_exp*q_exp(jm))
       else
         chii_cul=3.5D-3*(ti_m(jm)*1.D3*rlti_cul*tau_cul*falpha_cul)/
     &          (dsqrt(amassgas_exp)*bt_exp*q_exp(jm)*kappa_d**3)
       endif
c
       chietem=chie1_cul+chie2_cul
       chietim=0.D0
       chienem=0.D0
       chiitem=0.D0
c
       chiitim=chii_cul
       chiinem=0.D0
       difftem=0.D0
       difftim=0.D0
       diffnem=0.D0 
c  
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)  
c
c      write(*,100) jm, rho(jm), chietem, chiitim
c100   format(i2,2x,0p1f4.2,1p8e13.5)
c
      endif
c
      if (imodel.eq.7) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cjek Mixed-Shear model (Bohm-gyrobohm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Generalization of model proposed by Taroni, et al
c    in Plasma Phys. Control. Fusion 36 (1994) 1628.
c      chi_e = chi_e^B + chi_e^gB
c      chi_i = chi_i^B + chi_i^gB + chi_i^neo
c    Note: model include dependence on magnetic shear as proposed
c    by Romanelli and Zonca (Phys. Fluids B, 5 (1993) 4081)
c    Units: T (eV), B (T), radius (m), chi(m2/s)
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
       rlte_ms  = - zpmte
       rlti_ms  = - zpmti
       gfac=1.D0
       if(igeo_m.ge.1) gfac=geofac(jm)

       chietem=cmodel_e*gfac
     &       *abs( facbr_ms*rlte_ms*dbohm_ms*q_exp(jm)**2.D0
     &       *fshea_ms + facgb_ms*rlte_ms*dbohm_ms*gyro_ms/arho_exp )
       chietim=0.D0
       chienem=0.D0
       chiitem=0.D0
       chiitim=cmodel_i*gfac
     &       *abs( ceqei_ms*facbr_ms*rlti_ms*dbohm_ms*q_exp(jm)**2.D0
     &       *fshea_ms + facgb_ms*rlti_ms*dbohm_ms*gyro_ms/arho_exp )
       chiinem=0.D0
       difftem=0.D0
       difftim=0.D0
       diffnem=0.D0
c
       chietem=dmax1(chietem,chiemi_ms)
c
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
c...diagnostic printout
c
        if (lprintin_kb.ge.1) then
        open(4,file='mixedshear.out',status='unknown',
     &      access='sequential')
        write (4,*)
     &   'Variables used in mixed shear model'
        write (4,*)
        write (4,*) 'Rlte = ',rlte_ms
        write (4,*) 'Rlti = ',rlti_ms
        write (4,*) 'q = ',q_exp(jm)
        write (4,*) 'shat = ',shat_exp(jm)
        write (4,*) 'fshea = ',fshea_ms
        write (4,*) 'dbohm = ',dbohm_ms
        write (4,*) 'gyro = ',gyro_ms
        write (4,*) 'gfac = ',gfac
        write (4,*) 'chie = ',chietem
        write (4,*) 'chii = ',chiitem
        endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((imodel.eq.6).OR.(imodel.eq.66)) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Multi-mode 1995 model
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c JEK begin 9/15/98
c
c     Model of Weiland and Nordman
c     The old subroutine called is etawn8 from Bateman.
c     The new electromagnetic version w/ effects of finite beta and
c     parallel ion motion (hydrogenic and impurity) is etaw17a
c     Diagnostic ouput to etaw17a.out if open statement (nout_wn=1)
c     in etaw16.basf used and if lprint set larger in mlt0in.
c
c...  Note that:
c...  epsilon_n_wn = Ln/Lb as in the NF paper to compare
c...  the stability condition (Fig. 1) is zpmnh_wn,
c...  assuming that Lb = R**2/R0, given B=B0R0/R.
c     (JAK 3/23/95)
c
c     Changes:
c...  JEK 9/15/98 trying etaw17a, added wexb shearing rate to argument list
c...  JEK 9/15/98 subsituted in etaw16 instead of etaw14
c         added elongation and cetain_wn(25) to argument list
c...  JEK 12/02/96 replaced etaw12 by etaw14 (imodel=6) and switched
c     to normalized gradients
c
c...  JAK 3/30/95 elongation added to gradient lengths
c
c...  JAK 4/07/95 impurity effects included.
c     Note thatc zpmnh is not used in the main routine,
c     so it calculated here
c
c...  JAK 7/17/95 etawn8 replaced by etaw12
c...  Temporary: keep etaw8 (imodel=66) to alow consistency check on
c       etaw12
c
c
c       Model options
c          neq            effects
c           6      Hydrogenic + impurity ITG + TEM
c           7      6 eqns + collisions
c           8      6 eqns + collisions + hydrogenic parallel ion motion
c           9      6 eqns + collisions + finite beta
c                 + parallel ion motion in strong ballooning limit
c          10      6 eqns + collisions + finite beta
c                 + hydrogenic parallel ion motion
c          11      6 eqns + collisions + finite beta
c                 + hydrogenic and impurity parallel ion motion
c
c       Switches for Weiland advanced fluid model
c
c       cetain_wn(10) - coeff of k_\parallel in etaw17a
c       cetain_wn(11) - normalization for A_\parallel
c       cetain_wn(15) - coeff of nuhat in etaw17a
c       cetain_wn(20) - coeff of finite beta in etaw17a
c       cetain_wn(25) - coeff for elongation factor in Alfven frequency
c       cetain_wn(27) - coeff for ExB shearing rate
c       cetain_wn(30) - finite difference used to construct matrix
c       cetain_wn(32) - tolerance in eigenvalue solver
c       letain_wn(6) - set to zero for NAG eigenvalue solver
c       letain_wn(7) - equal zero to compute convective velocities
c       letain_wn(10) - set to zero to force complex
c       letain_wn(29) - greater than zero to print frequencies, fluxes
c
 888    letain_wn(6)=letain_wn6
        letain_wn(7)=letain_wn7
c       letain_wn(9)=2 ! use effective diffusivities
        cetain_wn(10)=1.D0
        cetain_wn(11)=1.D0
        cetain_wn(15)=cetain_wn15
        cetain_wn(20)=cetain_wn20
        cetain_wn(25)=cetain_wn25
        cetain_wn(27)=cetain_wn27
        cetain_wn(30)=cetain_wn30
        cetain_wn(32)=cetain_wn32
        nout_wn=1
c
        chietem=0.D0
        chiitem=0.D0
        diffnem=0.D0
c
c...constants
        zckb_kb = 1.60210D-16
        zce_kb  = 1.60210D-19
        zcmp_kb = 1.67252D-27
c
        fnsin_wn = nfast_exp(jm)/nem ! fast ion fraction ns/ne (JEK)
c
crew        zpmnh_wn = 
crew     >            -(dlog(ni_m(jm))-dlog(ni_m(jm+jd_m)))
crew     >            /(rho(jm)-rho(jm+jd_m))
        zpmnh_wn=zpmni
c
c       zpmti=8.D0/drhodr(jm)
c
        sqrtelong_wn = dsqrt(elong_exp(jm))
        if(igeo_m.ge.1) sqrtelong_wn=drhodr(jm)
        zpmns_wn = 1.D+3
        zpmne_wn = zpmne
        zpmte_wn = zpmte
        zpmth_wn = zpmti
        zpmtz_wn = zpmti
c
c...    calculate nz either from Zeff or from czin
c
        zeff_wn = zeff_exp(jm)
        if (zeff_wn.GE.1.0) then
          if (ABS(nem-nim).LT.1.0D-6) nim=nem-1.0D-6
c         czin_wn = (nem*zeff_wn-nim)/(nem-nim)
          czin_wn = 6.D0   !  carbon
          azin_wn = 12.D0 / amassgas_exp
c         nz_wn = (nem - nim)/czin_wn
          nz_wn = nz_exp(jm)
        else
          nz_wn = (nem - nim)/czin_wn
          zeff_wn = nim + nz_wn*(czin_wn**2.D0)/nem
        end if
        tz_wn = tim
c
        zpmnz_wn = 1.0D6
        if (nz_wn.GT.1.0D-6) then
         zpmnz_wn = (zpmne*nem - zpmnh_wn*nim)/(nem-nim)
        endif
c
        epsnein_wn =arho_exp/(sqrtelong_wn*rmaj_exp(jm)*zpmne_wn)   !L_ne/R
        epsnhin_wn =arho_exp/(sqrtelong_wn*rmaj_exp(jm)*zpmnh_wn)   !L_nH/R
        epsnzin_wn =arho_exp/(sqrtelong_wn*rmaj_exp(jm)*zpmnz_wn)   !L_nZ/R
        epsnsin_wn =arho_exp/(sqrtelong_wn*rmaj_exp(jm)*zpmns_wn)   !L_ns/R
        epstein_wn =arho_exp/(sqrtelong_wn*rmaj_exp(jm)*zpmte_wn)   !L_Te/R
        epsthin_wn =arho_exp/(sqrtelong_wn*rmaj_exp(jm)*zpmth_wn)   !L_Th/R
        epstzin_wn =arho_exp/(sqrtelong_wn*rmaj_exp(jm)*zpmtz_wn)   !L_Tz/R
c
c       epsnzin_wn = (nz_wn / nem) * czin_wn * epsnhin_wn
c    &   / ( (nz_wn / nem) * czin_wn - 1.D0 + epsnhin_wn / epsnein_wn )
c
c...normalized gradients JEK 12/02/96
c sign(max(abs(grdne(jz)),zepslon),grdne(jz)+zepsqrt)
        gnein_wn=1.D0/sign(dmax1(abs(epsnein_wn),1.D-6),epsnein_wn)   !R/L_ne
        gnhin_wn=1.D0/sign(dmax1(abs(epsnhin_wn),1.D-6),epsnhin_wn)   !R/L_nH
        gnzin_wn = 1.D0/sign(dmax1(abs(epsnzin_wn),1.D-6),epsnzin_wn) !R/L_nZ
        gnsin_wn = 1.D0/sign(dmax1(abs(epsnsin_wn),1.D-6),epsnsin_wn) !R/L_ns
        gtein_wn = 1.D0/sign(dmax1(abs(epstein_wn),1.D-6),epstein_wn) !R/L_Te
        gthin_wn = 1.D0/sign(dmax1(abs(epsthin_wn),1.D-6),epsthin_wn) !R/L_Th
        gtzin_wn = 1.D0/sign(dmax1(abs(epstzin_wn),1.D-6),epstzin_wn) !R/L_Tz
c
c       write(*,50) jm, rho(jm), gthin_wn, zpmth_wn
c
c...diagnostic arrays
c
        gne_wn(jm)=gnein_wn
        gnh_wn(jm)=gnhin_wn	
        gnz_wn(jm)=gnzin_wn	
        gns_wn(jm)=gnsin_wn
        gte_wn(jm)=gtein_wn
        gth_wn(jm)=gthin_wn
        gtz_wn(jm)=gtzin_wn
c
        tauhin_wn   = tim / tem                       ! Th/Te
        tauzin_wn   = tz_wn / tem                     ! Tz/Te
        fnzin_wn    = nz_wn / nem                     ! nZ/ne
c
cjak    czin_wn  = Z impurity charge number
cjak    azin_wn  = mZ/mH impurity mass to hydrogen isotope mass
c
        betaein_wn = 402.63*nem*tem/(1.D5*bt_exp**2.D0)  ! neTe/(B^2/2mu0)
        betahin_wn = 402.63*nim*tim/(1.D5*bt_exp**2.D0)  ! nHTH/(B^2/2mu0)
        betazin_wn = 402.63*nz_wn*tz_wn/(1.D5*bt_exp**2.D0)! nZTZ/(B^2/2mu0)
c
cjek    ftrapein_wn = 0.5  ! fraction of trapped electrons
        ftrapein_wn = dsqrt((2.D0*rho(jm)*arho_exp/rmaj_exp(jm))/
     >    (1.D0+rho(jm)*arho_exp/rmaj_exp(jm)))
c
c...    calculate rho_s as done in derivedexprofiles
c...    units are [cm]; convert to [m]
c...    Note that Z=1 is assumed in this formula.
c
        zgyrfi_wn  = zce_kb * bt_exp / (zcmp_kb*amassgas_exp)
        cs_wn   = dsqrt(zckb_kb*tem/(zcmp_kb*amassgas_exp))
        rhos_wn = cs_wn / zgyrfi_wn
        ky_wn = ekyrhoin_wn/rhos_wn
c
        omega_de_wn = (2.D0*ky_wn*rhos_wn*cs_wn)/rmaj_exp(jm)
c
c...    electron effective collision frequency (NLR standard)
c...    nueff = vu_th/(r/R)/w_de
c
        lnlam=24.D0-dlog(dsqrt(1.D14*nem)/1000.D0/tem)
        lnlam=dmax1(lnlam,1.D0)
        nue=2.91D-6*1.D14*nem*lnlam/(1000.D0*tem)**1.5D0
        vefin_wn = nue / (rho(jm)*arho_exp/rmaj_exp(jm)) / omega_de_wn
c
c       vefin_wn=1.D10  ! adiabatic electron test only
c
c...    q-profile, shear, and elongation
c
        qin_wn   = q_exp(jm)
        shearin_w= shat_exp(jm)
         if(igeo_m.ge.2.D0) shearin_w= shat_exp(jm)*drhodrrrho(jm)
        shearin_w = dmax1 ( 0.5D0, shearin_w )
        kappain_w= elong_exp(jm)
c
        ekparlin_wn = 1.D0 / dmax1 ( 1.D-6, qin_wn*gnhin_wn )
c
c...    ExB shearing rate
c       wexb_wn=0.D00
        wexb_wn=cetain_wn27*egamma_exp(jm)*csda_exp(jm)/omega_de_wn
        egamma_wn(jm) = wexb_wn
c
        if (imodel.EQ.66) then
         fig(1)=1.D0
         fig(2)=1.D0
         fig(3)=1.D0
        endif !imodel=66 
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c JEK Weiland model w/ up to 11 eqns using etaw17a
c
         if (imodel.eq.6) then
c
c         letain_wn(2) = 4
c         letain_wn(10)= 0
          fig(1)=0.8D0
          fig(2)=0.8D0
          fig(3)=0.8D0
         endif !imodel=6
         if(imodel.eq.6.or.imodel.eq.66) then
         if (lprintin_wn.ge.1)
     > open(1,file='etaw17a.out',status='unknown',access='sequential')
c
c... variables indicated by % are new compared to etawn8
c    Note: weiland14 used in NTCC module version of MM95
c          set itest_wn > 0. 
c
        if (itest_wn .gt. 0) then
c
        pi_m = datan2(0.D0, -1.D0)
        zckb_kb  = 1.60210D-16
        zcmu0_kb = 4.0D-7*pi_m
        zepslon_kb = 1.D-34
        zgasa_mm=apgasa
        czin_wn=apimpz
        azin_wn=apimpa / zgasa_mm
        nh_mm=ni_m(jm)*1.D19
        nz_mm=nz_m(jm)*1.D19
        nitot_mm=nh_mm+nz_mm
        zaimass_mm=( zgasa_mm * nh_mm + apimpa * nz_mm ) /
     &             nitot_mm
        btor_mm   = rmaj_exp(0)*bt_exp/rmaj_exp(jm)
        zep_kb    = dmax1( rmin_exp(jm)/rmaj_exp(jm), zepslon_kb)
        vthi_kb   = dsqrt(2.D0*zckb_kb*tim/(zcmp_kb*zaimass_mm))
        zgyrfi_kb = zce_kb*btor_mm/(zcmp_kb*zaimass_mm)
        zlari_kb  = vthi_kb/zgyrfi_kb
        zlarpo_kb = dmax1(zlari_kb*q_exp(jm)/zep_kb, zepslon_kb)
        zgmax_kb  = rmaj_exp(jm) / zlarpo_kb
        cs_wn   = dsqrt(zckb_kb*tem/(zcmp_kb*zaimass_mm))
        rhos_wn = cs_wn / zgyrfi_kb
        ky_wn = ekyrhoin_wn/rhos_wn
        omega_de_wn = (2.D0*ky_wn*rhos_wn*cs_wn)/rmaj_exp(jm)
        ftrapein_wn = dsqrt (2.D0*rmin_exp(jm)/
     &                (rmaj_exp(jm)*(1.D0+rmin_exp(jm)/rmaj_exp(jm) )))
        betaein_wn = 2.D0*zcmu0_kb*zckb_kb*nem*1.D19*tem/btor_mm**2
c
        drho = rho(jm-1)-rho(jm)
        zpmne_wn = -(dlog(ne_exp(jm-1))-dlog(ne_exp(jm)))/ drho
        zpmnh_wn = -(dlog(ni_exp(jm-1))-dlog(ni_exp(jm)))/ drho
        zpmnz_wn = -(dlog(nz_exp(jm-1))-dlog(nz_exp(jm)))/ drho
        zpmns_wn = -(dlog(nfast_exp(jm-1))-dlog(nfast_exp(jm)))/ drho
        zpmni_wn = -(dlog(nitot_exp(jm-1))-dlog(nitot_exp(jm)))/drho
        zpmnq_wn =  (dlog(q_exp(jm-1))-dlog(q_exp(jm)))/ drho
        zpmte_wn = -(dlog(te_exp(jm-1))-dlog(te_exp(jm)))/ drho
        zpmth_wn = -(dlog(ti_exp(jm-1))-dlog(ti_exp(jm)))/ drho
c        zpmth_wn = 4.0
        zpmtz_wn = zpmth_wn
c
        gnein_wn = drhodr(jm)*rmaj_exp(jm)*zpmne_wn/arho_exp ! R/L_ne
        gnhin_wn = drhodr(jm)*rmaj_exp(jm)*zpmnh_wn/arho_exp ! R/L_nH
        gnzin_wn = drhodr(jm)*rmaj_exp(jm)*zpmnz_wn/arho_exp ! R/L_nZ
        gnsin_wn = drhodr(jm)*rmaj_exp(jm)*zpmns_wn/arho_exp ! R/L_ns
        gniin_wn = drhodr(jm)*rmaj_exp(jm)*zpmni_wn/arho_exp ! R/L_ni
        gnqin_wn = drhodr(jm)*rmaj_exp(jm)*zpmnq_wn/arho_exp ! R/L_q
        gthin_wn = drhodr(jm)*rmaj_exp(jm)*zpmth_wn/arho_exp ! R/L_Th
        gtein_wn = drhodr(jm)*rmaj_exp(jm)*zpmte_wn/arho_exp ! R/L_Te
c
        gnein_wn =sign(dmin1(dabs(gnein_wn),zgmax_kb),gnein_wn)
        gnhin_wn =sign(dmin1(dabs(gnhin_wn),zgmax_kb),gnhin_wn)
        gnzin_wn =sign(dmin1(dabs(gnzin_wn),zgmax_kb),gnzin_wn)
        gnsin_wn =sign(dmin1(dabs(gnsin_wn),zgmax_kb),gnsin_wn)
        gniin_wn =sign(dmin1(dabs(gniin_wn),zgmax_kb),gniin_wn)
        gnqin_wn =sign(dmin1(dabs(gnqin_wn),zgmax_kb),gnqin_wn)
        gthin_wn =sign(dmin1(dabs(gthin_wn),zgmax_kb),gthin_wn)
        gtein_wn =sign(dmin1(dabs(gtein_wn),zgmax_kb),gtein_wn)
        gtzin_wn = gthin_wn
c
        shear_wn = gnqin_wn * rmin_exp(jm)/rmaj_exp(jm)
        zscyl_wn = dmax1(shear_wn,1.D-34)
        shearin_w = dmax1(0.5D0,zscyl_wn)
c
c       write(*,50) jm, rho(jm), gthin_wn, zpmth_wn
c       write(*,50) jm, rho(jm), czin_wn, azin_wn,
c    &              apimpa, zaimass_mm
c
         call weiland14 (
     &   letain_wn,     !i* input and control variable
                        !  letain(2) number of elements computed for
                        !            transport matrix
                        !  letain(6) > 1 use IMSL to compute eigenvalues
                        !            <=1 use NAG to compute eigenvalues
                        !  letain(7) = 0 compute convective velocities
                        !            alternatively, set the convective
                        !            velocities to zero and rescale
                        !            the diffusivity matrix
                        !  letain(10) = 0 force complex valued
                        !            eigenvalue matrices
     &   cetain_wn,     !i* input and control variable
                        !  cetain(10) coeff. of kpar (default to 1.0)
                        !  cetain(30) finite difference used to construct
                        !         transport matrix (see zgm(j1,jd) matrix)
                        !  cetain(32) tolerance used in eigenvalue solver
     &   lprintin_wn,   !i* controls printout (higher values more printout)
                        !  lprintin_wn=-1 means no output at all
     &   neq_wn,        !i* number of equations
     &   nout_wn,       !i% output lun
     &   gnein_wn,      !i R/L_ne
     &   gnhin_wn,      !i R/L_nH
     &   gnzin_wn,      !i R/L_nZ
     &   gtein_wn,      !i R/L_Te
     &   gthin_wn,      !i R/L_Th
     &   gtzin_wn,      !i R/L_Tz
     &   tauhin_wn,     !i Th/Te
     &   tauzin_wn,     !i Tz/Te
     &   fnzin_wn,      !i nZ/ne
     &   czin_wn,       !i/o* Z impurity charge number -> zimpz
     &   azin_wn,       !i/o* mZ/mH impurity mass to hydrogen isotope mass
     &   fnsin_wn,      !i ns/ne fraction of superthermal hydrogenic ions
     &   betaein_wn,    !i% neTe/(B^2/2mu0)
     &   ftrapein_wn,   !i fraction of trapped electrons
     &   vefin_wn,      !i% nueff = nuth/(r wDe/R)
     &   qin_wn,        !i% q-value
     &   shearin_w,     !i% shear = d ln q / d ln r
     &   ekyrhoin_wn,   !i* normalized poloidal wave number
     &   wexb_wn,       !i* ExB shearing rate
     &   ndim_wn,       !i* first dimension of the 2-D array difthi
                        !  and the maximum number of unstable modes allowed
     &   omega_wn,      !o real*8 part of the frequencies
                        !  normalized by $ \omega_{De} $
     &   gamma_wn,      !o growth rates (normalized)
     &   difthi_wn,     !  diffusivity matrix
                        !  normalized by $ k_y^2 / \omega_{De} $
     &   velthi_wn,     !o convective velocities
                        !  normalized by $ k_y / \omega_{De} $
     &   chieff_wn,     !o effective total diffusivities for nhTh,
                        !  nh, neTe, nZ, nZTZ
                        !  normalized by $ k_y^2 / \omega_{De} $
     &   flux_wn,       !o flux
     &   nmodes_wn,     !o number of unstable modes
     &   nerr)          !o%
c
        else
c
        call etaw17a (
     &   letain_wn,     !i* input and control variable
                        !  letain(2) number of elements computed for
                        !            transport matrix
                        !  letain(6) > 1 use IMSL to compute eigenvalues
                        !            <=1 use NAG to compute eigenvalues
                        !  letain(7) = 0 compute convective velocities
                        !            alternatively, set the convective
                        !            velocities to zero and rescale
                        !            the diffusivity matrix
                        !  letain(10) = 0 force complex valued
                        !            eigenvalue matrices
     &   cetain_wn,     !i* input and control variable
                        !  cetain(10) coeff. of kpar (default to 1.0)
                        !  cetain(30) finite difference used to construct
                        !         transport matrix (see zgm(j1,jd) matrix)
                        !  cetain(32) tolerance used in eigenvalue solver
     &   lprintin_wn,   !i* controls printout (higher values more printout)
                        !  lprintin_wn=-1 means no output at all
     &   neq_wn,        !i* number of equations
     &   nout_wn,       !i% output lun
     &   gnein_wn,      !i R/L_ne
     &   gnhin_wn,      !i R/L_nH
     &   gnzin_wn,      !i R/L_nZ
     &   gtein_wn,      !i R/L_Te
     &   gthin_wn,      !i R/L_Th
     &   gtzin_wn,      !i R/L_Tz
     &   tauhin_wn,     !i Th/Te
     &   tauzin_wn,     !i Tz/Te
     &   fnzin_wn,      !i nZ/ne
     &   czin_wn,       !i/o* Z impurity charge number
     &   azin_wn,       !i/o* mZ/mH impurity mass to hydrogen isotope mass
     &   fnsin_wn,      !i ns/ne fraction of superthermal hydrogenic ions
     &   betaein_wn,    !i% neTe/(B^2/2mu0)
     &   betahin_wn,    !i nHTH/(B^2/2mu0)
     &   betazin_wn,    !i/o nZTZ/(B^2/2mu0)
     &   ftrapein_wn,   !i fraction of trapped electrons
     &   vefin_wn,      !i% nueff = nuth/(r wDe/R)
     &   qin_wn,        !i% q-value
     &   shearin_w,     !i% shear = d ln q / d ln r
     &   kappain_w,     !i% elongation
     &   ekyrhoin_wn,   !i* normalized poloidal wave number
     &   ekparlin_wn,   !i* kparr Ln
     &   wexb_wn,       !i* ExB shearing rate
     &   ndim_wn,       !i* first dimension of the 2-D array difthi
                        !  and the maximum number of unstable modes allowed
     &   omega_wn,      !o real*8 part of the frequencies
                        !  normalized by $ \omega_{De} $
     &   gamma_wn,      !o growth rates (normalized)
     &   difthi_wn,     !  diffusivity matrix
                        !  normalized by $ k_y^2 / \omega_{De} $
     &   velthi_wn,     !o convective velocities
                        !  normalized by $ k_y / \omega_{De} $
     &   chieff_wn,     !o effective total diffusivities for nhTh,
                        !  nh, neTe, nZ, nZTZ
                        !  normalized by $ k_y^2 / \omega_{De} $
     &   nmodes_wn,     !o number of unstable modes
     &   perform_wn,    !o performance index from IMSL gpirg
     &   nerr)          !o%
c
         endif          !itest_wn=0 or otherwise
         endif          !imodel=6 or imodel=66
c
c...   Note that in the definitions of Bateman:
c...
c...      d( (n~ +T~)v e )/dt = chieff_wn * gradTe
c...      where n~ is the perturbed n, etc.
c...
c...      But the heat flux is 3/2 d(nT)/dt = 3/2 d( (n~+T~)v~e )/dt
c...
c...      so chie = 3/2 * chieff_wn
c...
c...       nd his actual power balance equation reads:
c...
c...      d(nT)/dt = 2/3(nabla . q_e) + 2/3 sourceterms
c
c      gfac=1.
c      if(igeo_m.ge.1) gfac=geofac(jm)
c      chietem = gfac*chieff_wn(3) * omega_de_wn/(ky_wn**2.D0) * 1.5D0
c      chietim = 0.D0
c      chienem = 0.D0
c      chiitem = 0.D0
c      chiitim = gfac*chieff_wn(1) * omega_de_wn/(ky_wn**2.D0) * 1.5D0
c      chiinem = 0.D0
c      difftem = 0.D0
c      difftim = 0.D0
c      diffnem = gfac*chieff_wn(2) * omega_de_wn/(ky_wn**2.D0) * 1.D0
c
c      if (incl_elong_kb.eq.1) then    ! --- elongation scaling ---
c      elong_kb = elong_exp(jm)
c       if (elong_kb.ne.0) then
c        chietem = chietem*(elong_kb**(pow_elong_kb))
c        chietim = chietim*(elong_kb**(pow_elong_kb))
c        chienem = chienem*(elong_kb**(pow_elong_kb))
c        chiitem = chiitem*(elong_kb**(pow_elong_kb))
c        chiitim = chiitim*(elong_kb**(pow_elong_kb))
c        chiinem = chiinem*(elong_kb**(pow_elong_kb))
c        difftem = difftem*(elong_kb**(pow_elong_kb))
c        difftim = difftim*(elong_kb**(pow_elong_kb))
c        diffnem = diffnem*(elong_kb**(pow_elong_kb))
c       endif
c      endif !(incl_elong_kb=1)
c
c      chiegb_m(jm)=chietem/cgyrobohm_m(jm)
c      chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
c      diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
c new stuff
       gfac=1.D0
       if(igeo_m.ge.6) gfac=geofac(jm)
       if(igeo_m.ge.8) gfac=gradrhosq_exp(jm)
       diff_wn_kb=fig(1)*gfac*chieff_wn(2)*
     >            omega_de_wn/(ky_wn**2.D0)*1.D0
       chie_wn_kb=fig(2)*gfac*chieff_wn(3)*
     >            omega_de_wn/(ky_wn**2.D0)*1.0D0
       chii_wn_kb=fig(3)*gfac*chieff_wn(1)*
     >            omega_de_wn/(ky_wn**2.D0)*1.0D0
c
c check on gamma ordering
c possible dimension error if nmaxwn > ikymax_gf
       gamma_test_max=-1.D6
       do ik=1,nmaxwn
        gamma_k_j(ik,jm)=gamma_wn(ik)*omega_de_wn/csda_m(jm)
        freq_k_j(ik,jm)=omega_wn(ik)*omega_de_wn/csda_m(jm)
        if(gamma_k_j(ik,jm).gt.gamma_test_max) then
          gamma_test_max=gamma_k_j(ik,jm)
          anrate_m(jm)=gamma_k_j(ik,jm) 
          anfreq_m(jm)=freq_k_j(ik,jm)
        endif
       enddo
c
c       write(*,'(i3,2x,0p1f4.2,1p8e12.4)') jm, rho(jm),
c    &       zpmti, zpmti*drhodr(jm), gthin_wn, chii_wn_kb
c       write(*,'(i3,2x,0p1f4.2,1p8e12.4)') jm, rho(jm),
c    &       anrate_m(jm), gamma_wn(3), 
c    &       gamma_wn(3)*omega_de_wn/csda_m(jm), 
c    &       omega_wn(3), omega_de_wn, csda_m(jm)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Add in Guzdar-Drake model
      if (incl_gd_kb.eq.1) then
c
c constants
c
        zckb_kb = 1.60210D-16
        zcme_kb = 9.1091D-31
        zcmp_kb = 1.67252D-27
        zce_kb  = 1.60210D-19
        zcc_kb = 2.997925D+8
        zcmu0_kb   = 4.0D-7 * pi_m
        zceps0_kb  = 1.D0 / ( zcc_kb**2.D0 * zcmu0_kb )
        zepslon_kb = 1.D-34
c
      if (itest_wn .gt. 0) then
        frb(1)=1.D0
        frb(2)=1.D0
        frb(3)=1.D0
        zgddia_kb = 0.15D0
        nh_mm=ni_m(jm)*1.D19
        nz_mm=nz_m(jm)*1.D19
        nitot_mm=nh_mm+nz_mm
        btor_mm    = rmaj_exp(0)*bt_exp/rmaj_exp(jm)
        zep_kb     = dmax1( rmin_exp(jm)/rmaj_exp(jm), zepslon_kb)
        vthe_kb    = dsqrt(2.D0*zckb_kb*tem/zcme_kb)
        vthi_kb    = dsqrt(2.D0*zckb_kb*tim/(zcmp_kb*zaimass_mm))
        zgyrfe_kb  = zce_kb*bt_exp/zcme_kb
        zgyrfi_kb  = zce_kb*btor_mm/(zcmp_kb*zaimass_mm)
        zlare_kb   = vthe_kb/zgyrfe_kb
        zlari_kb   = vthi_kb/zgyrfi_kb
        zlarpo_kb   = dmax1(zlari_kb*q_exp(jm)/zep_kb, zepslon_kb)
        zlog_kb    = 37.8D0-dlog(dsqrt(nem*1.D19)/tem)
        zcf_kb     = (4.D0*dsqrt(pi_m)/3.D0)*
     >               (zce_kb/(4.D0*pi_m*zceps0_kb))**2.D0
        zcf_kb     = zcf_kb*(zce_kb/zckb_kb)*
     >               dsqrt((zce_kb/zcme_kb)*(zce_kb/zckb_kb))
        znuei_kb   = zcf_kb*dsqrt(2.D0)*nem*1.D19*zlog_kb*
     >               zeff_exp(jm)/(tem*dsqrt(tem))
        zgmax_kb   = rmaj_exp(jm) / zlarpo_kb
c
c  zgpr = -R ( d p   / d r ) / p    for thermal pressure
c
        zgpr_kb = ( nem*1.D19 * tem * ( gnein_wn + gtein_wn )
     &            + nitot_mm  * tim * ( gniin_wn + gthin_wn ) )
     &            / ( nem*1.D19 * tem + nitot_mm * tim )
        zgpr_kb = sign ( dmin1 ( dabs ( zgpr_kb ), zgmax_kb ), zgpr_kb )
        zgdp_kb = 2.D0*pi_m*((q_exp(jm)*zlare_kb)**2.D0)*znuei_kb*
     &           zgpr_kb*100.D0*zgddia_kb
c
        gfac = 1.D0
        if(igeo_m.ge.6) gfac=geofac(jm)
        if(igeo_m.ge.8) gfac=gradrhosq_exp(jm)
c
        diff_gd_kb = frb(1) * gfac * zgdp_kb
        chie_gd_kb = frb(2) * gfac * zgdp_kb
        chii_gd_kb = frb(3) * gfac * zgdp_kb
c
      else
c
        frb(1)=15.D0
        frb(2)=15.D0
        frb(3)=15.D0
c
        vthe_kb    = dsqrt(2.D0*zckb_kb*tem/zcme_kb)
        zgyrfe_kb  = zce_kb * bt_exp / zcme_kb
        zlare_kb   = vthe_kb / zgyrfe_kb
        vthi_kb    = dsqrt(2.D0*zckb_kb*tim/(zcmp_kb*amassgas_exp))
        zgyrfi_kb  = zce_kb * bt_exp / (zcmp_kb*amassgas_exp)
        zlari_kb   = vthi_kb / zgyrfi_kb
        zlog_kb    = 37.8D0-dlog(dsqrt(nem*1.D19)/tem)
        zcf_kb     = (4.D0*dsqrt(pi_m)/3.D0)*(zce_kb/
     >               (4.D0*pi_m*zceps0_kb))**2.D0
     >               *(zce_kb/zckb_kb)*zce_kb/dsqrt(zcme_kb*zckb_kb)
        znuei_kb   = zcf_kb*dsqrt(2.D0)*nem*1.D19*zlog_kb*zeff_exp(jm)
     >              /(tem*dsqrt(tem))
c
c       inverse aspect ratio
        zep_kb     = dmax1( arho_exp*rho(jm)/rmaj_exp(jm),zepslon_kb)
c
c...    Note that zlpr is the max of a read in value, or zlarpo.
c...    The read in value is zslpr(jz) which is Lp [cm]
c...    Include sqrt(kappa).
c
        zlarpo_kb  = dmax1(zlari_kb*q_exp(jm)/zep_kb,zepslon_kb)
        zsglpr_kb  = sign ( 1.D0, zlpr_kb )
cjak    zlpr_kb  =arho_exp*(1.D0/zpmte + 1.D0/zpmne)*100. ![cm]
cjek    zlpre_kb =arho_exp*(1.D0/(sqrtelong_wn*zpmte) +1.D0/
cjek >    (sqrtelong_wn*zpmne))
cjek    zlpri_kb =arho_exp*(1.D0/(sqrtelong_wn*zpmti) + 1.D0/
cjek >    (sqrtelong_wn*zpmni))
cjek    zlpr_kb    = zlpre_kb + zlpri_kb ![m]
cjek    zlpr_kb    = dmax1(abs(zlpr_kb),zlarpo_kb)*zsglpr_kb
c
c...jek new pressure gradient scale length
c
        sqrtelong_kb = dsqrt(elong_exp(jm))
        if(igeo_m.ge.1) sqrtelong_kb=drhodr(jm)
c
        zepsne =arho_exp/(sqrtelong_kb*rmajor_exp*zpmne)   !L_ne/R
        zepsni =arho_exp/(sqrtelong_kb*rmajor_exp*zpmni)   !L_ni/R
        zepste =arho_exp/(sqrtelong_kb*rmajor_exp*zpmte)   !L_Te/R
        zepsti =arho_exp/(sqrtelong_kb*rmajor_exp*zpmti)   !L_Ti/R
c
        zgne =1.D0/sign(dmax1(dabs(zepsne),1.D-6),zepsne)  !R/L_ne
        zgni =1.D0/sign(dmax1(dabs(zepsni),1.D-6),zepsni)  !R/L_ni
        zgte =1.D0/sign(dmax1(dabs(zepste),1.D-6),zepste)  !R/L_Te
        zgti =1.D0/sign(dmax1(dabs(zepsti),1.D-6),zepsti)  !R/L_Th
c
        zprth   = nem*tem + nim*tim
        zdprth  = nem*tem*(zgne+zgte) + nim*tim*(zgni+zgti)
        zsgdpr  = sign ( 1.D0, zdprth)
        zdpdr   = dmax1(abs(zdprth),zlarpo_kb) * zsgdpr
        zlpr_kb = dmax1(rmajor_exp*zprth/zdpdr,1.D-6)
c
c...end new stuff
c
        zsound_kb  = dsqrt(zckb_kb*tem/(zcmp_kb*amassgas_exp))
        zgdtime_kb = dsqrt ( rmaj_exp(jm)*zlpr_kb/2.D0 )/zsound_kb
        zgdrlp_kb  = (rmaj_exp(jm)/zlpr_kb)**.25D0         ! (R/L_p)^.25
        zgdlna_kb  = 2.D0 * pi_m * q_exp(jm) * zgdrlp_kb
        zrhos_kb   = zsound_kb/zgyrfi_kb
        zgdln_kb   = zgdlna_kb*dsqrt((znuei_kb*rmaj_exp(jm)*zrhos_kb)
     >              /(2.D0*zgyrfe_kb)) ! char length
c
c       fixed diamagnetic stabilization (alpha set to zero, zgddia_kb = 1)
c
        zgdalf_kb = (zrhos_kb*zsound_kb*zgdtime_kb) / (zlpr_kb*zgdln_kb)
        zgddia_kb = 1.D0 / ( 1.D0 + zgdalf_kb**2.D0 )
        zgddia_kb = 1.D0
        zgdp_kb    = 2.D0 * pi_m * (q_exp(jm) * zlare_kb)**2.D0*znuei_kb
     >              * ( rmaj_exp(jm) / zlpr_kb ) * zgddia_kb
c
        zdtite_kb  = exp( -pow_tem_kb * ((tim/tem)-1.D0)**2.D0 )
c
        gfac=1.D0
        if(igeo_m.ge.6) gfac=geofac(jm)
        if(igeo_m.ge.8) gfac=gradrhosq_exp(jm)
        diff_gd_kb = frb(1) * gfac * zgdp_kb * zdtite_kb
        chie_gd_kb = frb(2) * gfac * zgdp_kb * zdtite_kb
        chii_gd_kb = frb(3) * gfac * zgdp_kb * zdtite_kb
c
c Diagnostic printout
c
        if (lprintin_gd.ge.1) then
        open(2,file='guzdrake.out',status='unknown',access='sequential')
        write (2,*)
     &   'Variables used in Guzdar-Drake model'
        write (2,*)
        write (2,*) 'amass = ',amassgas_exp
        write (2,*) 'R = ',rmaj_exp(jm)
        write (2,*) 'Bt = ',bt_exp
        write (2,*) 'Zeff = ',zeff_exp(jm)
        write (2,*) 'epsilon = ',zep_kb
        write (2,*) 'ne = ',nem
        write (2,*) 'ni = ',nim
        write (2,*) 'Te = ',tem
        write (2,*) 'Ti = ',tim
        write (2,*) 'vthe = ',vthe_kb
        write (2,*) 'vthi = ',vthi_kb
        write (2,*) 'cs = ',zsound_kb
        write (2,*) 'omegae = ',zgyrfe_kb
        write (2,*) 'omegai = ',zgyrfi_kb
        write (2,*) 'rhoe = ',zlare_kb
        write (2,*) 'rhoi = ',zlari_kb
        write (2,*) 'rhos = ',zrhos_kb
        write (2,*) 'zlarpo = ',zlarpo_kb
        write (2,*) 'loglambda = ',zlog_kb
        write (2,*) 'zcf = ',zcf_kb
        write (2,*) 'nuei = ',znuei_kb
        write (2,*) 'arho = ',arho_exp
        write (2,*) 'Lte = ',zepste*rmajor_exp
        write (2,*) 'Lne = ',zepsne*rmajor_exp
        write (2,*) 'Lpr = ',zlpr_kb
        write (2,*) 'q = ',q_exp(jm)
        write (2,*) 'ideal growth rate = ',zgdtime_kb
        write (2,*) 'zgdrlp = ',zgdrlp_kb
        write (2,*) 'zgdlna_kb = ',zgdlna_kb
        write (2,*) 'char length = ',zgdln_kb
        write (2,*) 'zgdp = ',zgdp_kb
        write (2,*) 'zdtite = ',zdtite_kb
        write (2,*) 'gfac = ',gfac
        write (2,*) 'chi_gd = ',chie_gd_kb
        endif
c
      endif
      else
        diff_gd_kb=0.D0
        chie_gd_kb=0.D0
        chii_gd_kb=0.D0
      endif
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Add in Singer-Rewoldt-Bateman kinetic ballooning model
      if (incl_sr_kb.eq.1) then
        fkb(1)=1.D0
        fkb(2)=0.65D0
        fkb(3)=0.65D0
c       cthery8_kb=3.5D0
        cthery3_kb=0.5D0
c
        pi_m = datan2(0.D0, -1.D0)
        zcmu0_kb   = 4.0D-7 * pi_m
        zckb_kb = 1.60210D-16
        zcmp_kb = 1.67252D-27
        zepslon_kb = 1.D-34
        zlgeps_kb = dlog ( zepslon_kb )
c
      if (itest_wn .gt. 1) then
        btor_mm    = rmaj_exp(0)*bt_exp/rmaj_exp(jm)
        zsound_kb  = dsqrt(zckb_kb*tem/(zcmp_kb*zaimass_mm))
        zgyrfi_kb  = zce_kb*btor_mm/(zcmp_kb*zaimass_mm)
        zrhos_kb   = zsound_kb/zgyrfi_kb
        nh_mm=ni_m(jm)*1.D19
        nz_mm=nz_m(jm)*1.D19
        nitot_mm=nh_mm+nz_mm
        zgpr_kb = ( nem*1.D19 * tem * ( gnein_wn + gtein_wn )
     &            + nitot_mm * tim * ( gniin_wn + gthin_wn ) )
     &            / ( nem*1.D19 * tem + nitot_mm * tim )
        zgpr_kb = sign(dmin1 ( dabs ( zgpr_kb ), zgmax_kb ), zgpr_kb)
        zbeta_kb  = (2.D0 * zcmu0_kb * zckb_kb / btor_mm**2) * 
     &              (nem*1.D19 * tem + nitot_mm * tim)
        zbprim_kb = dabs( zbeta_kb * zgpr_kb / rmaj_exp(jm) )
        qin_kb    = q_exp(jm)
        shear_wn  = gnqin_wn*rmin_exp(jm)/rmaj_exp(jm)
        zscyl_wn  = dmax1(shear_wn,1.D-34)
        shearin_w = dmax1(0.5D0,zscyl_wn)
        zbc1_kb   = dabs(shearin_w)/(1.7D0*qin_kb**2.D0*rmaj_exp(jm))
        zfbthn_kb = dexp( dmin1(dabs(zlgeps_kb),
     &     dmax1(-dabs(zlgeps_kb),3.5D0*(zbprim_kb/zbc1_kb - 1.D0))) )
        zdk_kb    = dabs( zsound_kb*zrhos_kb**2.D0* 
     &              zfbthn_kb*zgpr_kb/rmaj_exp(jm) )
c
        gfac=1.D0
        if(igeo_m.ge.6) gfac=geofac(jm)
        if(igeo_m.ge.8) gfac=gradrhosq_exp(jm)
c
        if ( zgpr_kb .gt. 0.D0 ) then
          diff_sr_kb = fkb(1) * gfac * zdk_kb
          chie_sr_kb = fkb(2) * gfac * zdk_kb
          chii_sr_kb = fkb(3) * gfac * zdk_kb
        else
          diff_sr_kb = 0.D0
          chie_sr_kb = 0.D0
          chii_sr_kb = 0.D0
        endif
c
      else
c
        zsound_kb  = dsqrt(zckb_kb*tem/(zcmp_kb*amassgas_exp))
        zgyrfi_kb  = zce_kb * bt_exp / (zcmp_kb*amassgas_exp)
        zrhos_kb   = zsound_kb/zgyrfi_kb
c
c..pressure gradient scale length Lp [m]
c
        sqrtelong_kb = dsqrt(elong_exp(jm))
        if(igeo_m.ge.1) sqrtelong_kb=drhodr(jm)
c
        zepsne =arho_exp/(sqrtelong_kb*rmaj_exp(jm)*zpmne)   !L_ne/R
        zepsni =arho_exp/(sqrtelong_kb*rmaj_exp(jm)*zpmni)   !L_ni/R
        zepste =arho_exp/(sqrtelong_kb*rmaj_exp(jm)*zpmte)   !L_Te/R
        zepsti =arho_exp/(sqrtelong_kb*rmaj_exp(jm)*zpmti)   !L_Ti/R
c
        zgne =1.D0/sign(dmax1(dabs(zepsne),1.D-6),zepsne)  !R/L_ne
        zgni =1.D0/sign(dmax1(dabs(zepsni),1.D-6),zepsni)  !R/L_ni
        zgte =1.D0/sign(dmax1(dabs(zepste),1.D-6),zepste)  !R/L_Te
        zgti =1.D0/sign(dmax1(dabs(zepsti),1.D-6),zepsti)  !R/L_Th
c
        zprth   = nem*tem + nim*tim
        zdprth  = nem*tem*(zgne+zgte) + nim*tim*(zgni+zgti)
        zsgdpr  = sign ( 1.D0, zdprth)
        zdpdr   = dmax1(dabs(zdprth),zlarpo_kb) * zsgdpr
        zlpr_kb = dmax1(rmaj_exp(jm)*zprth/zdpdr,1.D-4)
c
c..total beta
c
        betaein_kb = 402.63D0*nem*tem/(1.D5*bt_exp**2.D0)  ! neTe/(B^2/2mu0)
        betaiin_kb = 402.63D0*nim*tim/(1.D5*bt_exp**2.D0)  ! niTi/(B^2/2mu0)
        zbeta_kb=betaein_kb + betaiin_kb
c
c..q-profile and shear
c
        qin_kb   = q_exp(jm)
        shearin_k= dmax1(cthery3_kb,shat_exp(jm))
        if(igeo_m.ge.2) then
          shearin_k= dmax1(cthery3_kb,shat_exp(jm)*drhodrrrho(jm))
        endif
c
c..beta primes and crits
c
        zbprim_kb = dabs(zbeta_kb/zlpr_kb)
        zbcoef1_kb = 1.D0
        zbcoef2_kb = 1.D0
        zbc1_kb = zbcoef1_kb*dabs(shearin_k)/
     &            (1.7D0*qin_kb**2*rmaj_exp(jm))
        zbpbc1 = zbprim_kb /zbc1_kb
        zbc2_kb = zbcoef2_kb*4.D0*dabs(shearin_k)/
     &            (qin_kb**2*rmaj_exp(jm))
c
        zfbthn_kb = dexp( dmin1(dabs(zlgeps_kb),
     &     dmax1(-dabs(zlgeps_kb),cthery8_kb*(zbprim_kb/zbc1_kb-1.D0))))
        zdk_kb = dabs( zsound_kb * zrhos_kb**2 * zfbthn_kb / zlpr_kb )
c
        zdtite_kb  = dexp( -pow_tem_kb * ((tim/tem)-1.D0)**2.D0 )
c
        gfac=1.D0
        if(igeo_m.ge.6) gfac=geofac(jm)
        if(igeo_m.ge.8) gfac=gradrhosq_exp(jm)
        diff_sr_kb = fkb(1) * gfac * zdk_kb * zdtite_kb
        chie_sr_kb = fkb(2) * gfac * zdk_kb * zdtite_kb
        chii_sr_kb = fkb(3) * gfac * zdk_kb * zdtite_kb
c
c       write(*,31) jm, rho(jm), chii_sr_kb*elong_exp(jm)**pow_elong_kb
c       write(*,51) jm, rho(jm), zlare_kb, znuei_kb, zlpr_kb
 31    format(i2,2x,0p1f9.6,1p2e14.6,' KB transport')
c
c from theory.f:
c     zbetae = 2. * zcmu0 * zckb * zne * zte /   zb**2
c     zbeta=(2.*zcmu0*zckb/zb**2.D0)*(zne*zte+zni*zti)
c     zbprim = abs(zbeta/zlpr)
c     zbcoef1 = 1.0
c     zbcoef2 = 1.0
c     if ( abs(cthery(78)) .gt. zepslon ) zbcoef1 = cthery(78)
c     if ( abs(cthery(79)) .gt. zepslon ) zbcoef2 = cthery(79)
c     zbc1   = zbcoef1 * abs(zshat)/(1.7*zq**2*zrmaj)
c     zbpbc1 = zbprim/zbc1
c     zbc2   = zbcoef2 * 4.* abs(zshat)/(zq**2*zrmaj)
c     zelfkb = zelonf**cthery(15)
c       zfbthn = exp( min(abs(zlgeps),
c    &     dmax1(-dabs(zlgeps),cthery(8)*(zbprim/zbc1 - 1.D0))) )
c       zdk = dabs( zsound * zrhos**2 * zfbthn / zlpr )
c       zdkb(jz)  = zdk*fkb(1)*zelfkb
c       thkbe(jz) = zdk*fkb(2)*zelfkb
c       thkbi(jz) = zdk*fkb(3)*zelfkb
c
        if (lprintin_kb.ge.1) then
        open(4,file='kinball.out',status='unknown',access='sequential')
        write (4,*)
     &   'Variables used in kinetic ballooning model'
        write (4,*)
        write (4,*) 'amass = ',amassgas_exp
        write (4,*) 'R = ',rmaj_exp(jm)
        write (4,*) 'Bt = ',bt_exp
        write (4,*) 'ne = ',nem
        write (4,*) 'ni = ',nim
        write (4,*) 'Te = ',tem
        write (4,*) 'Ti = ',tim
        write (4,*) 'Lpr = ',zlpr_kb
        write (4,*) 'cs = ',zsound_kb
        write (4,*) 'rhos = ',zrhos_kb
        write (4,*) 'beta = ',zbeta_kb
        write (4,*) 'q = ',qin_kb
        write (4,*) 'shear = ',shearin_k
        write (4,*) 'beta prime = ',zbprim_kb
        write (4,*) 'beta prime c1 = ',zbc1_kb
        write (4,*) 'beta prime c2 = ',zbc2_kb
        write (4,*) 'zlgeps = ',zlgeps_kb
        write (4,*) 'zfbthn = ',zfbthn_kb
        write (4,*) 'zdk = ',zdk_kb
        write (4,*) 'zdtite = ',zdtite_kb
        write (4,*) 'gfac = ',gfac
        endif
c
      endif
      else
        diff_sr_kb=0.D0
        chie_sr_kb=0.D0
        chii_sr_kb=0.D0
      endif
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Add up total transport
       chietem = chie_wn_kb + chie_gd_kb + chie_sr_kb
       chietim = 0.D0
       chienem = 0.D0
       chiitem = 0.D0
       chiitim = chii_wn_kb + chii_gd_kb + chii_sr_kb
       chiinem = 0.D0
       difftem = 0.D0
       difftim = 0.D0
       diffnem = diff_wn_kb + diff_gd_kb + diff_sr_kb
c
       if (incl_elong_kb.eq.1) then    ! --- elongation scaling ---
       elong_kb = elong_exp(jm)
        if (elong_kb.ne.0) then
         chietem = chietem*(elong_kb**(pow_elong_kb))
         chietim = chietim*(elong_kb**(pow_elong_kb))
         chienem = chienem*(elong_kb**(pow_elong_kb))
         chiitem = chiitem*(elong_kb**(pow_elong_kb))
         chiitim = chiitim*(elong_kb**(pow_elong_kb))
         chiinem = chiinem*(elong_kb**(pow_elong_kb))
         difftem = difftem*(elong_kb**(pow_elong_kb))
         difftim = difftim*(elong_kb**(pow_elong_kb))
         diffnem = diffnem*(elong_kb**(pow_elong_kb))
        endif
       endif !(incl_elong_kb=1)

       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c      if (i_proc.eq.0) then
c        write(*,50) jm, rho(jm), chii_wn_kb*elong_kb**pow_elong_kb,
c    &             chii_sr_kb*elong_kb**pow_elong_kb,
c    &             chii_gd_kb*elong_kb**pow_elong_kb,
c    &             chiitim,cgyrobohm_m(jm),zpmth_wn
c      endif
c
      endif
c
c end Weiland or Multi-mode model
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (imodel.eq.67) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Multi-mode 1995 model (NTCC version)
c 11-01/00 J. Kinsey
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c...constants
        zckb_kb = 1.60210D-16
        zce_kb  = 1.60210D-19
        zcmp_kb = 1.67252D-27
c
        rmin_mm=rmin_exp(jm)
        rmaj_mm=rmaj_exp(jm)
        elong_mm=elong_exp(jm)
        ne_mm=ne_m(jm)*1.D19
        nh_mm=ni_m(jm)*1.D19
        nz_mm=nz_m(jm)*1.D19
        ns_mm=ns_m(jm)*1.D19
        nitot_mm=nh_mm+nz_mm
        nse_mm=0.D0  ! electron density from fast ions
        te_mm=te_m(jm)
        ti_mm=ti_m(jm)
        q_mm=q_exp(jm)
        btor_mm=rmaj_exp(0)*bt_exp/rmaj_exp(jm)
        zimpz_mm=apimpz
        zimpa_mm=apimpa
        zgasz_mm=apgasz
        if (idata.eq.0 .and. (nm1_exp(jm)+nm2_exp(jm)).gt.1.D-3) then
          zgasa_mm=( 1.D0*nm1_exp(jm) + apgasa*nm2_exp(jm) )/
     &             ( nm1_exp(jm)  + nm2_exp(jm) )
        else
          zgasa_mm=apgasa
        endif
        zaimass_mm=( zgasa_mm * nh_mm + zimpa_mm * nz_mm ) /
     &             nitot_mm
c       zeff_mm=zeff_exp(jm)
        zeff_mm=( nh_mm * zgasz_mm**2.D0 + nz_mm*zimpz_mm**2.D0 +
     &            ns_mm * abgasz**2.D0) / ne_mm   
        sqrtelong_wn = dsqrt(elong_exp(jm))
c
c... Diamagnetic frequency
        zgyrfi_wn = zce_kb * bt_exp / (zcmp_kb*amassgas_exp)
        cs_wn = dsqrt(zckb_kb*tem/(zcmp_kb*amassgas_exp))
        rhos_wn = cs_wn / zgyrfi_wn
        ky_wn = ekyrhoin_wn/rhos_wn
        omega_de_wn = (2.D0*ky_wn*rhos_wn*cs_wn)/(rmaj_exp(jm))
c
c... ExB shear rate
        wexb_wn=cetain_wn27*egamma_exp(jm)*csda_exp(jm)/omega_de_wn
        egamma_wn(jm) = wexb_wn
c
c... Calculate normalized gradients:
c
c  zgrdne(jz) = - R ( d n_e / d r ) / n_e
c  zgrdnh(jz) = - R ( d n_h / d r ) / n_h
c     n_h = thermal hydrogenic density (sum over hydrogenic species)
c  zgrdnz(jz) = - R ( d Z n_Z / d r ) / ( Z n_Z )
c     n_Z = thermal impurity density,  Z = average impurity charge
c           summed over all impurities
c  zgrdns(jz) = - R ( d Z_s n_s / d r ) / ( Z_s n_s )
c     n_s = superthermal ion density,
c     Z_s = average superthermal ion charge at gridpoint jz
c  zgrdni(jz) = - R ( d n_i / d r ) / n_i
c     n_i = thermal ion density (sum over hydrogenic and impurity)
c  zgrdte(jz) = - R ( d T_e / d r ) / T_e
c  zgrdti(jz) = - R ( d T_i / d r ) / T_i
c  zgrdq (jz) =   R ( d q   / d r ) / q    related to magnetic shear
c
c  Note: nz calculated from Zeff in rep_iter.f if not provided
c        in ufiles. also, tz set equal to ti
c
        drho = rho(jm-1)-rho(jm)
        zpmne_wn = -(dlog(ne_exp(jm-1))-dlog(ne_exp(jm)))/ drho
        zpmnh_wn = -(dlog(ni_exp(jm-1))-dlog(ni_exp(jm)))/ drho
        zpmnz_wn = -(dlog(nz_exp(jm-1))-dlog(nz_exp(jm)))/ drho
        zpmns_wn = -(dlog(nfast_exp(jm-1))-dlog(nfast_exp(jm)))/ drho
        zpmni_wn = -(dlog(nitot_exp(jm-1))-dlog(nitot_exp(jm)))/drho
        zpmnq_wn =  (dlog(q_exp(jm-1))-dlog(q_exp(jm)))/ drho
        zpmte_wn = -(dlog(te_exp(jm-1))-dlog(te_exp(jm)))/ drho
        zpmth_wn = -(dlog(ti_exp(jm-1))-dlog(ti_exp(jm)))/ drho
        zpmtz_wn = zpmth_wn
c
        gnein_wn = drhodr(jm)*rmaj_exp(jm)*zpmne_wn/arho_exp ! R/L_ne
        gnhin_wn = drhodr(jm)*rmaj_exp(jm)*zpmnh_wn/arho_exp ! R/L_nH
        gnzin_wn = drhodr(jm)*rmaj_exp(jm)*zpmnz_wn/arho_exp ! R/L_nZ
        gnsin_wn = drhodr(jm)*rmaj_exp(jm)*zpmns_wn/arho_exp ! R/L_ns
        gniin_wn = drhodr(jm)*rmaj_exp(jm)*zpmni_wn/arho_exp ! R/L_ni
        gnqin_wn = drhodr(jm)*rmaj_exp(jm)*zpmnq_wn/arho_exp ! R/L_q
        gthin_wn = drhodr(jm)*rmaj_exp(jm)*zpmth_wn/arho_exp ! R/L_Th
        gtein_wn = drhodr(jm)*rmaj_exp(jm)*zpmte_wn/arho_exp ! R/L_Te
c
        gnein_wn =sign(dmax1(dabs(gnein_wn),1.D-6),gnein_wn)
        gnhin_wn =sign(dmax1(dabs(gnhin_wn),1.D-6),gnhin_wn)
        gnzin_wn =sign(dmax1(dabs(gnzin_wn),1.D-6),gnzin_wn)
        gnsin_wn =sign(dmax1(dabs(gnsin_wn),1.D-6),gnsin_wn)
        gniin_wn =sign(dmax1(dabs(gniin_wn),1.D-6),gniin_wn)
        gnqin_wn =sign(dmax1(dabs(gnqin_wn),1.D-6),gnqin_wn)
        gthin_wn =sign(dmax1(dabs(gthin_wn),1.D-6),gthin_wn)
        gtein_wn =sign(dmax1(dabs(gtein_wn),1.D-6),gtein_wn)
        gtzin_wn = gthin_wn
c
c...diagnostic arrays
c
        gne_wn(jm)=gnein_wn
        gnh_wn(jm)=gnhin_wn	
        gnz_wn(jm)=gnzin_wn	
        gns_wn(jm)=gnsin_wn
        gte_wn(jm)=gtein_wn
        gth_wn(jm)=gthin_wn
        gtz_wn(jm)=gtzin_wn
c
        ndim_wn=5    ! matrix dimension
        npoints=1    ! number of radial pts
        nnout=44     ! unit number for diagnostic printout
        lprint_wn=0  ! no diagnostic printout
        lsuper=0     ! 0=std, 1 for supershot settings
        lreset=0     ! use default internal settings for model
c
c...initialize switches
c
        do jr=1,25
          cswitch(jr) = 0.D0
        enddo
        do jr=1,8
          lswitch(jr) = 0
        enddo
c
        call mmm95(
     &   rmin_mm,  rmaj_mm,  elong_mm
     & , ne_mm,    nh_mm,    nz_mm,    nse_mm
     & , zeff_mm,  te_mm,    ti_mm,    q_mm,       btor_mm
     & , zimpz_mm, zimpa_mm, zgasa_mm, zaimass_mm, wexb_wn
     & , gnein_wn, gniin_wn, gnhin_wn, gnzin_wn
     & , gtein_wn, gthin_wn, gnqin_wn
     & , zthiig,   zthdig,   ztheig,   zthzig
     & , zthirb,   zthdrb,   ztherb,   zthzrb
     & , zthikb,   zthdkb,   zthekb,   zthzkb
     & , zgamma_wn,zomega_wn,difthi_wn,velthi_wn,  flux_wn
     & , ndim_wn,  npoints,  nnout,    lprint_wn,  nerr
     & , lsuper,   lreset,   lswitch,  cswitch
     & , fig_mm,   frb_mm,   fkb_mm )
c
       chietem = ztheig + ztherb + zthekb
       chiitim = zthiig + zthirb + zthikb
       diffnem = zthdig + zthdrb + zthdkb
c
       gfac=1.D0
c      if(igeo_m.ge.1) gfac=geofac(jm)
       chietem = chietem*gfac
       chiitim = chiitim*gfac
       diffnem = diffnem*gfac
c
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
c      write(*,51) jm, rho(jm), zthiig, zthirb,
c    &             zthikb
c
      endif
c
c end MMM95 model
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (imodel.eq.8) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 2DGLF quasilinear model  GLF23
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
c  test effect of canonical density gradients
c note pzmn_sol only in zpmne_q and zpmni_q which go into
c glf2d drivers rlne_gf and rlni_gf

       zpmne_q=(1.D0+xparam_gf(20))*zpmne*pzmn_sol(jm)-xparam_gf(19)*
     > (dlog(q_exp(jm-1))-dlog(q_exp(jm)))/(rho(jm-1)-rho(jm))
       zpmni_q=(1.D0+xparam_gf(20))*zpmni*pzmn_sol(jm)-xparam_gf(19)*
     > (dlog(q_exp(jm-1))-dlog(q_exp(jm)))/(rho(jm-1)-rho(jm))


       rmaj_gf=rmaj_exp(jm)/arho_exp
       rmin_gf=rmin_exp(jm)/arho_exp
       if(igeo_m.eq.-1) then
        rmaj_gf=rmajor_exp/arho_exp
        rmin_gf=rho(jm)*arho_exp/sqrt(elonga_exp)
       endif
       q_gf=q_exp(jm)
       betae_gf=cbetae*betae_m(jm)
       shat_gf=shat_exp(jm)
        if(igeo_m.ge.2) shat_gf=shat_exp(jm)*drhodrrrho(jm)
       alpha_gf=xalpha*alpha_exp(jm)

      if(ialphastab.eq.2) then
       if(igeo_m.ge.0) then
c2/27/98        alpha_m(jm)=sqrt(elong_exp(jm))/
c               ((1.D0+elong_exp(jm)**2.D0)/2.D0)*
         alpha_m(jm)=drhodr(j)*
     >    q_exp(jm)**2.D0*rmaj_exp(jm)/arho_exp*
     >    betae_m(jm)*((ti_m(jm)/te_m(jm)*ni_m(jm)/ne_m(jm))*
     >   (zpmni+zpmti)
     >    +zpmne+zpmte)
c2/27/98     >    betae_m(jm)*((ti_m(jm)/te_m(jm))*(zpmne+zpmti)
c2/27/98     >    +zpmne+zpmte)  
       endif
       if(igeo_m.ge.4) then
        alpha_m(jm)=geoalpha(j)*
     >    q_exp(jm)**2.D0*rmaj_exp(jm)/arho_exp*
     >    betae_m(jm)*((ti_m(jm)/te_m(jm)*ni_m(jm)/ne_m(jm))*
     >   (zpmni+zpmti)
     >    +zpmne+zpmte)
c2/27/98     >    betae_m(jm)*((ti_m(jm)/te_m(jm))*(zpmne+zpmti)
c2/27/98     >    +zpmne+zpmte) 
       endif
       if(igeo_m.eq.-1) then
        alpha_m(jm)=
     > q_exp(jm)**2.D0*rmajor_exp/arho_exp*
     > sqrt(elonga_exp)/((1.D0+elonga_exp**2.D0)/2.D0)*
     >    betae_m(jm)*((ti_m(jm)/te_m(jm)*ni_m(jm)/ne_m(jm))*
     >   (zpmni+zpmti)
     >    +zpmne+zpmte)
c2/27/98     >    betae_m(jm)*((ti_m(jm)/te_m(jm))*(zpmne+zpmti)
c2/27/98     >    +zpmne+zpmte) 
       endif
      endif   
       if(ialphastab.gt.0) then
        alpha_gf=xalpha*alpha_m(jm)
       endif
       elong_gf=elong_exp(jm)
       if(igeo_m.eq.-1) elong_gf=elonga_exp
       xnu_gf=cxnu*xnu_m(jm)
       taui_gf=tim/tem
       amassgas_gf=amassgas_exp

       rlte_gf=zpmte*sqrt(elong_exp(jm))
       rlti_gf=zpmti*sqrt(elong_exp(jm))
       rlne_gf=zpmne_q*sqrt(elong_exp(jm))
       rlni_gf=rlne_gf
       if(i_dengrad.ge.1) then
        rlne_gf=zpmne_q*sqrt(elong_exp(jm))
        rlni_gf=zpmni_q*sqrt(elong_exp(jm)) 
       endif
       dil_gf=0.D0
       apwt_gf=1.D0
       aiwt_gf=0.D0
       rlnimp_gf=1.D0
       zpmnimp=1.D0
       if(i_dengrad.eq.2) dil_gf=1.D0-nim/nem
       if(i_dengrad.eq.3) then
        apwt_gf=nim/nem
        aiwt_jp1=(zeff_exp(jm+1)*ne_exp(jm+1)-ni_exp(jm+1)
     >    -nfast_exp(jm+1))/(zimp_gf**2.D0*ne_exp(jm+1))
        xnimp_jp1=aiwt_jp1*ne_exp(jm+1)
        aiwt_gf=(zeff_exp(jm)*ne_exp(jm)-ni_exp(jm)
     >    -nfast_exp(jm))/(zimp_gf**2.D0*ne_exp(jm))
        xnimp=aiwt_gf*ne_exp(jm)
        
        zpmnimp=-(dlog(xnimp_jp1)-dlog(xnimp))/(rho(jm+1)-rho(jm))
        rlnimp_gf=zpmnimp*sqrt(elong_exp(jm))
c zimp_gf and amassimp_gf are inputs
       endif


        if(igeo_m.ge.1) then
         rlte_gf=zpmte*drhodr(jm)
         rlti_gf=zpmti*drhodr(jm)
         rlne_gf=zpmne_q*drhodr(jm)
         rlni_gf=rlne_gf
         rlnimp_gf=zpmnimp*drhodr(jm)

         if(i_dengrad.ge.1) then
          rlne_gf=zpmne_q*drhodr(jm)
          rlni_gf=zpmni_q*drhodr(jm)
         endif
        endif
        if(igeo_m.eq.-1) then
         rlte_gf=zpmte*sqrt(elonga_exp)
         rlti_gf=zpmti*sqrt(elonga_exp)
         rlnimp_gf=zpmnimp*sqrt(elonga_exp)
         rlne_gf=zpmne_q*sqrt(elonga_exp)
         rlni_gf=rlne_gf
         if(i_dengrad.ge.1) then
          rlne_gf=zpmne_q*sqrt(elonga_exp)
          rlni_gf=zpmni_q*sqrt(elonga_exp)
          rlnimp_gf=zpmnimp*sqrt(elonga_exp)
         endif
        endif


       gamma_star_gf=vstarp_exp(jm)
       gamma_e_gf=egamma_exp(jm)
       gamma_p_gf=gamma_p_exp(jm)
       gamma_mode_gf=gamma_mode_exp(jm)
       if(irotstab.gt.0) then
        gamma_star_gf=vstarp_m(jm)
        gamma_e_gf=egamma_m(jm)
        if(i_delay.gt.1) gamma_e_gf=egamma_d(jm,1)
        if(i_gridsmooth.gt.0) gamma_e_gf=egamma_g(jm)
        gamma_p_gf=gamma_p_m(jm)
        gamma_mode_gf=gamma_mode_m(jm)
       endif

      if(irotstab.eq.2) then
       vstar_sign=-1.D0
        j=jm
        vstarp_m(jm)=
     > +vstar_sign*(rho(j+1)+rho(j))/2.D0*(
     >(ti_m(j)/te_m(j))*csda_m(j)*(zpni_m(j+1))
     >   *rhosda_m(j)/rho(j+1)
     >-(ti_m(j)/te_m(j))*csda_m(j)*(zpni_m(j))
     >   *rhosda_m(j)/rho(j)
     >  )/(rho(j+1)-rho(j))/csda_m(j)
     >  -vstar_sign*zpmti*(ti_m(j)/te_m(j))
     >  *csda_m(j)*(zpni_m(j))*rhosda_m(j)/csda_m(j) 
  
        gamma_star_gf=vstarp_m(jm)  

      endif
  
        call glf2d
  
c 4/25/96 note: chie_gf and chii_gf from glf2d are energy diffusivies
c but here we can redefine them as heat diffusivites if we
c take the convection from experiment using xconv=1 or 5./3.

       chie_gf=chie_gf-xconv*1.5D0*rlni_gf/rlte_gf*diff_gf
       chii_gf=chii_gf-xconv*1.5D0*rlni_gf/rlti_gf*diff_gf

       ky_j(jm)=xkyf_gf
       gamma_j(jm,1)=gamma_gf(1)
       gamma_j(jm,2)=gamma_gf(2)
       gamma_j(jm,3)=gamma_gf(3)
       gamma_j(jm,4)=gamma_gf(4)
  
       freq_j(jm,1)=freq_gf(1)
       freq_j(jm,2)=freq_gf(2)
       freq_j(jm,3)=freq_gf(3)
       freq_j(jm,4)=freq_gf(4)
 
       phi_norm_j(jm,1)=phi_norm_gf(1)
       phi_norm_j(jm,2)=phi_norm_gf(2)
       phi_norm_j(jm,3)=phi_norm_gf(3)
       phi_norm_j(jm,4)=phi_norm_gf(4)
  
       do ik=1,ikymax_gf
        gamma_k_j(ik,jm)=gamma_k_gf(1,ik)
        freq_k_j(ik,jm) =freq_k_gf(1,ik)
        chie_k_j(ik,jm) = chie_k_gf(ik)
        chii_k_j(ik,jm) = chii_k_gf(ik)
       enddo

       anrate_m(jm)=gamma_gf(1)
       dnrate_m(jm)=0.D0
       dtnrate_m(jm)=0.D0
       anfreq_m(jm)=freq_gf(1)
       dnfreq_m(jm)=0.D0
        gfac=1.
        if(igeo_m.ge.1) gfac=geofac(jm)
        if(igeo_m.eq.-1) 
     >     gfac=(1.D0+elonga_exp**2.D0)/
     >     (2.D0*elonga_exp)/gradrhosq_exp(jm)
       chiegb_m(jm)=chie_gf*cmodel*gfac
       chiigb_m(jm)=chii_gf*cmodel*gfac
       diffgb_m(jm)=diff_gf*cmodel*gfac
       diffgb_im_m(jm)=diff_im_gf*cmodel*gfac
c note diff_im_gf used only for diagnostics
c  to be compared with diff_gf or chii_gf

       exchgb_m(jm)=exch_gf*cmodel*gfac
       etagb_phi_m(jm)=eta_phi_gf*cmodel*gfac
     
       etagb_par_m(jm)=eta_par_gf*cmodel*gfac
       etagb_per_m(jm)=eta_per_gf*cmodel*gfac
    
       chie_e_gb_m(jm)=chie_e_gf*cmodel*gfac
  
c exch_m in MW/m**3  

      exch_m(jm)=1.D19*kevdsecpmw*nem*tem*csda_m(jm)*
     > rhosda_m(jm)**2.D0*exch_gf*cmodel

 
c exch_m is directly related to the flow
c for a single mode branch exch_gf=-(-freq_gf(1)/xkyf_gf)*diff_gf*rln_gf.
c we can  not expect to find exch_m without knowing flow_exp as input.
c and solving self consistant flow eq. flown=flow_exp for density
c density profile.
c
c however, knowing freq_gf(1) from the gf model we can compute exch_exp  
c from flow_exp using
c       flowm=kevdsecpmw*1.*nem*1.D19/arho_exp*gradrhosq_exp(jm)*
c     >       sfactor(jm)*(difftem*zpmte+difftim*zpmti+diffnem*zpmne)
c we have:

       diffgb_local=flow_exp(jm)/
     > (kevdsecpmw*1.D0*nem*1.D19/arho_exp*gradrhosq_exp(jm)*
     > sfactor(jm)*zpmne_q)/cgyrobohm_m(jm)
 
       exchgb_local=-(-freq_gf(1)/xkyf_gf)*diffgb_local*rlni_gf
  
       exch_exp(jm)=1.D19*
     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2.D0*exchgb_local
    
c       exch_exp(jm)=flow_exp(jm)*tem*(-1.D0)*(-freq_gf(1)/xkyf_gf)*
c     > sqrt(elong_exp(jm))*arho_exp/gradrhosq_exp(jm)/sfactor(jm)
c     >/arho_exp(jm)**2.D0
  

c   note electron(ion) wave freq > 0(<0) cool(heat) electrons
c (-1) denotes electron to ion
 
c to emphasize, we can not know exch_exp better than we know flow_exp



c chietem, chiitim, diffen in m**2/sec

      chietem=cmodel*gfac*chie_gf*cgyrobohm_m(jm)
      chietim=0.D0
      chienem=0.D0
 
      chiitem=0.D0
      chiitim=cmodel*gfac*chii_gf*cgyrobohm_m(jm)
      chiinem=0.D0
  
      difftem=0.D0
      difftim=0.D0
      diffnem=cmodel*gfac*diff_gf*cgyrobohm_m(jm)
   
   
      etaparm=cmodel*gfac*eta_par_gf*cgyrobohm_m(jm)
      etaperm=cmodel*gfac*eta_per_gf*cgyrobohm_m(jm)
  
      etaphim=cmodel*gfac*eta_phi_gf*cgyrobohm_m(jm)

      diff_m(jm)=diffnem
      chie_m(jm)=chietem
      chii_m(jm)=chiitim
      etapar_m(jm)=etaparm
      etaper_m(jm)=etaperm

      etaphi_m(jm)=etaphim
c
      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if (imodel.eq.81) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     imodel.eq.8 with self contained precall
c
      if(i_proc .eq. 0 .and. iparam_pt(20).eq.1) 
     >   write(6,*) 'before callglf2d'
c     iparam_pt7=iparam_pt(7)
c     iparam_pt8=iparam_pt(8)
      alpha_e=alpha_e_gf
      x_alpha=xalpha
c
      jshoot=0
c     jeigen=0    ! cgg solver
      igrad=1     ! input gradients
      if (iglf2d.gt.0) igrad=0
      if (iglf2d.lt.0) igrad=0
      if(iparam_pt(6).eq.-1) jshoot=0
      jmm=jm
c
c store input gradients when using DV method
c
      if (igrad .eq. 1) then
        zpmte_glf(jm)=zpmte
        zpmti_glf(jm)=zpmti
        zpmne_glf(jm)=zpmne
        zpmni_glf(jm)=zpmni
      endif
c
      do ik=1,ikymax_gf
       gamma_k_j(ik,jm)   = 0.D0
       freq_k_j(ik,jm)    = 0.D0
       diff_k_j(ik,jm)    = 0.D0
       chie_k_j(ik,jm)    = 0.D0
       chie_etg_j(ik,jm)  = 0.D0
       chie_itg_j(ik,jm)  = 0.D0
       chii_k_j(ik,jm)    = 0.D0
       etaphi_k_j(ik,jm)  = 0.D0
      enddo
c
c... onetwo testing
c
      if (iglf2d.lt.-1) jmm=0
c      zpmti=4.5
c
c     if calling during an intermediate step i_delay should be 0 or less
c
      call callglf2d(
     >                 !INPUTS
     > jeigen,         ! eigenvalue solver
     >                 ! 0 for cgg (default), 1 for tomsqz, 2 for zgeev
     > nroot_gf,       ! no. roots,8 for default, 12 for impurity dynamics
     > jshoot,         ! jshoot=0 time-dep code;jshoot=1 shooting code
     > jmm,            ! grid number;jmm=0 does full grid jm=1 to jmaxm-1
     > jmaxm,          ! profile grids 0 to jmaxm
     > itport_pt,      ! 1:5 transport flags
     > irotstab,       ! 0 to use egamma_exp; 1 use egamma_m
     > ave_ve,         ! smoothing parameter for ExB shear velocity (> 0 more smoothing)
     > te_m,           ! 0:jmaxm te Kev           itport_pt(2)=1 transport
     > ti_m,           ! 0:jmaxm ti Kev           itport_pt(3)=1 transport
     > ne_m,           ! 0:jmaxm ne 10**19 1/m**3
     > ni_m,           ! 0:jmaxm ni 10**19 1/m**3 itport_pt(1)=1 transport
     > ns_m,           ! 0:jmaxm ns 10**19 1/m**3 
     > igrad,          ! default 0, for D-V method use i_grad=1 to input gradients 
     > i_dengrad,      ! default 2, for simple dilution
     > zpmte,          ! override jm  log gradient te  iparam_pt8=-1
     > zpmti,          ! override jm  log gradient ti  iparam_pt8=-1
     > zpmne,          ! override jm  log gradient ne  iparam_pt8=-1
     > zpmni,          ! override jm  log gradient ni  iparam_pt8=-1
     > angrotp_exp,    ! 0:jmaxm exp plasma toroidal angular velocity 1/sec
     >                 ! if itport_pt(4)=0 itport_pt(5)=0
     > egamma_exp,     ! 0:jmaxm exp exb shear rate in units of csda_exp
     >                 ! if itport_pt(4)=-1 itport_pt(5)=0
     > gamma_p_exp,    ! 0:jmaxm exp par. vel. shear rate in units of csda_exp
     >                 ! if itport_pt(4)=-1 itport_pt(5)=0
     > vphi_m*corot,   ! 0:jmaxm toroidal velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=0 otherwise output
     > vpar_m,         ! 0:jmaxm parallel velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
     > vper_m,         ! 0:jmaxm perp. velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
     > zeff_exp,       ! 0:jmaxm ne in 10**19 1/m**3
     > bt_exp,         ! vaccuum axis toroidal field in tesla
     > bt_flag,        ! switch for effective toroidal field use in rhosda
     > bteff_exp,      ! 0:jmaxm effective toroidal field in tesla
     > rho,            ! 0:jmaxm 0 < rho < 1 toroidal flux surface label
     > arho_exp,       ! unit length toroidal flux surface LCFS in meters
     >                 !   toroidal flux= B0*rho_phys**2/2,
     >                 !   B0=bt_exp, arho_exp=rho_phys_LCFS
     > gradrho_exp,    ! 0:jmaxm dimensionless <|grad rho_phys |**2>
     > gradrhosq_exp,  ! 0:jmaxm dimensionless <|grad rho_phys |>
     >                 !NOTE:can set arho_exp=1.,if gradrho_exp=<|grad rho |>
     >                 !                 and gradrhosq_exp = <|grad rho |**2>
     > rmin_exp,       ! 0:jmaxm minor radius in meters
     > rmaj_exp,       ! 0:jmaxm major radius in meters
     > sfactor,        ! 0:jmaxm shape factor needed for gfac
     > vprime,         ! 0:jmaxm dvol/dr needed for gfac
     > rmajor_exp,     ! axis major radius
     > zimp_exp,       ! effective Z of impurity
     > amassimp_exp,   ! effective A of impurity
     > q_exp,          ! 0:jmaxm saftey factor
     > shat_exp,       ! 0:jmaxm  rho d q_exp/ d rho
     > alpha_exp,      ! 0:jmaxm MHD alpha from experiment
     > pfast_exp,      ! 0:jmaxm fast ion pressure
     > elong_exp,      ! 0:jmaxm  elongation
     > amassgas_exp,   !  atomic number working hydrogen gas
     > alpha_dia,      ! switch for diamagnetic term in ExB
     > alpha_e,        ! 1 full (0 no) no ExB shear stab
     > x_alpha,        ! 1 full (0 no) alpha stabilization  with alpha_exp
     >                 !-1 full (0 no) self consistent alpha_m stab.
     > i_delay,        !i_delay time delay for ExB shear should be non-zero only
     >                 ! once per step and is less than or equal 10
     >                 !OUTPUTS
     > diffnem,        ! ion plasma diffusivity in m**2/sec
     > chietem,        ! electron ENERGY diffuivity in m**2/sec
     > chiitim,        ! ion      ENERGY diffuivity in m**2/sec
     > etaphim,        ! toroidal velocity diffusivity in m**2/sec
     > etaparm,        ! parallel velocity diffusivity in m**2/sec
     > etaperm,        ! perpendicular velocity diffusivity in m**2/sec
     > exchm,          ! turbulent electron to ion ENERGY exchange in MW/m**3
     >                 ! 0:jmaxm values
     >  diff_m,
     >  chie_m,
     >  chii_m,
     >  etaphi_m,
     >  etapar_m,
     >  etaper_m,
     >  exch_m,
     >
     >  egamma_m,      !0:jmaxm exb shear rate in units of local csda_m
     >  egamma_d,      !0:jmaxm exb shear rate delayed by i_delay steps
     >  gamma_p_m,     !0:jmaxm par. vel. shear rate in units of local csda_m
     >  anrate_m,      !0:jmaxm leading mode rate in unints of local csda_m
     >  anrate2_m,     !0:jmaxm 2nd mode rate in units of local csda_m
     >  anfreq_m,      !0:jmaxm leading mode frequency
     >  anfreq2_m,     !0:jmaxm 2nd mode frequency
     >  xkymax_m,      !0:jmaxm leading mode ky
     >  xkymax2_m      !0:jmaxm 2nd leading mode ky
     > )    
       if(iparam_pt(20).eq.1) write(6,*) ' after callglf2d'
c
c... store contributions to transport from individual k's
c
       if (istep_glf.eq.lprint_glfstep .and. idvflag.eq.1) then
        do ik=1,ikymax_gf
         gamma_k_j(ik,jm)= gamma_k_gf(1,ik)
         freq_k_j(ik,jm) = freq_k_gf(1,ik)
         diff_k_j(ik,jm) = diff_k_gf(ik)*cgyrobohm_m(jm)*geofac(jm)
         chie_itg_j(ik,jm) = chie_k_gf(ik)*cgyrobohm_m(jm)*geofac(jm)     ! ITG/TEM
         chie_etg_j(ik,jm) = chie_e_k_gf(ik)*cgyrobohm_m(jm)*geofac(jm)   ! ETG
         chie_k_j(ik,jm) = chie_itg_j(ik,jm) +  chie_etg_j(ik,jm)         ! ITG/TEM & ETG
         chii_k_j(ik,jm) = chii_k_gf(ik)*cgyrobohm_m(jm)*geofac(jm)       ! ITG/TEM
         etaphi_k_j(ik,jm) = eta_phi_k_gf(ik)*cgyrobohm_m(jm)*geofac(jm)  ! ITG/TEM
c        write(*,'(i3,2x,i3,1p8e11.3)') jm, ik, xkyf_k_gf(ik), chietem,
c    >       chie_itg_j(ik,jm), chie_etg_j(ik,jm), chie_k_j(ik,jm)
c        write(*,'(i3,2x,i3,1p10e11.3)') jm, ik, xkyf_k_gf(ik),
c    >       gamma_k_j(ik,jm), freq_k_j(ik,jm),
c    >       chie_itg_j(ik,jm), chie_etg_j(ik,jm), chii_k_j(ik,jm)
c    >       diff_k_j(ik,jm), chie_k_j(ik,jm), chii_k_j(ik,jm), 
c    >       zpmne, zpmni
        enddo
c       write(*,'(i3,2x,1p10e11.3)') jm, anrate_m(jm), anfreq_m(jm), 
c    >                               zpmni, Psour(jm), diffnem
       endif
c
       chiegb_m(jm)=chietem/cgyrobohm_m(jm)
       chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
       diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
c
       etagb_phi_m(jm)=etaphi_m(jm)/cgyrobohm_m(jm)
       etagb_par_m(jm)=etapar_m(jm)/cgyrobohm_m(jm)
       etagb_per_m(jm)=etaper_m(jm)/cgyrobohm_m(jm)
c
       chiegb_glf(jm)=chietem
       chiigb_glf(jm)=chiitim
       diffgb_glf(jm)=diffnem
       etagb_phi_glf(jm)=etaphi_m(jm)
       etagb_par_glf(jm)=etapar_m(jm)
       etagb_per_glf(jm)=etaper_m(jm)
c
c     write(*,60) jm, rho(jm), anrate_m(jm), chie_m(jm), chii_m(jm)
c
c...new stuff
       if (jm.eq.lprint_glfgrid .and. istep_glf.eq.lprint_glfstep 
     >     .and. lprint_glf2.gt.0 .and. idvflag.eq.0) then
        open (35,file='glf2.out',status='unknown')
        write(*,53) lprint_glfgrid, rho(jm), lprint_glfstep
        write(*,55)
        do j=1,10
          write(*,60) j, xkyf_k_gf(j), gamma_k_gf(1,j),
     >                freq_k_gf(1,j), phi_norm_k_gf(1,j),
     >                chie_k_gf(j), chie_e_k_gf(j), chii_k_gf(j)
        enddo
        write(*,65) anrate_m(jm), anfreq_m(jm)
        write(*,70) chiegb_glf(jm),chiigb_glf(jm),etagb_phi_glf(jm)
        write(*,75) chie_gf,chie_e_gf,chii_gf
        write(*,*) zpmne_glf(jm)
        write(*,*) zpmni_glf(jm)
        write(*,*) q_exp(jm)
        write(*,*) shat_exp(jm)
        write(*,*) alpha_exp(jm)
        write(*,*) gamma_p_m(jm)
        write(*,*) chiitim
        close(35)
       endif
       idvflag=idvflag+1
       if(idvflag.eq.2) idvflag=0
c
c      if(idvflag.eq.0) then
c        write(*,50) jm, rho(jm), egamma_exp(jm), csda_exp(jm)
c      endif
c      write(*,50) jm, rho(jm), vphi_m(jm), bteff_exp(jm),
c    >             alpha_exp(jm), gamma_p_exp(jm)
c
c... printout for ONETWO testing
c
      if (istep_glf.eq.0.and.idvflag.eq.0.and.iglf2d.eq.99) then
      if (jm.eq.1) then
        write(*,'(a12,a9)') 'shot     = ',shot
        write(*,'(a12,i3)') 'leigen   = ',jeigen
        write(*,'(a12,i3)') 'iglf     = ',iglf
        write(*,'(a12,i3)') 'jmaxm    = ',jmaxm
        write(*,'(a12,5(x,i3))')'itport_pt = ',(itport_pt(i),i=1,5)
        write(*,'(a12,i3)') 'nroot    = ',nroot_gf
        write(*,'(a12,i3)') 'jshoot   = ',jshoot
        write(*,'(a12,i3)') 'jmm      = ',jmm
        write(*,'(a12,i3)') 'irotstab = ',irotstab
        write(*,'(a12,i3)') 'igrad    = ',igrad
        write(*,'(a12,i3)') 'idengrad = ',i_dengrad
        write(*,'(a12,i3)') 'bt_flag  = ',bt_flag
        write(*,'(a12,f10.6)') 'bt_exp   = ',bt_exp
        write(*,'(a12,f10.6)') 'arho_exp = ',arho_exp
        write(*,'(a12,f10.6)') 'amassgas = ',amassgas_exp
        write(*,'(a12,f10.6)') 'amassimp = ',amassimp_exp
        write(*,'(a12,f10.6)') 'zimp_exp = ',zimp_exp
        write(*,'(a12,f10.6)') 'alpha_e  = ',alpha_e
        write(*,'(a12,f10.6)') 'x_alpha  = ',x_alpha
      endif
      write(*,50) jm, rho(jm), te_m(jm), ti_m(jm), ne_m(jm), ni_m(jm),
     >            zeff_exp(jm), ns_m(jm)
      write(*,50) jm, rho(jm), rmin_exp(jm),
     >            rmaj_exp(jm), gradrho_exp(jm),
     >            gradrhosq_exp(jm), elong_exp(jm)
      write(*,50) jm, rho(jm), elong_exp(jm), q_exp(jm), shat_exp(jm),
     >            alpha_exp(jm), angrotp_exp(jm), egamma_exp(jm)
      write(*,50) jm, rho(jm), zpmne, zpmni, zpmte, zpmti,
     >            vphi_m(jm), egamma_m(jm)
      write(*,50) jm, rho(jm), diffnem, chietem, chiitim, exchm
      write(*,50) jm, rho(jm), etaphim, etaparm, etaperm
      write(*,50) jm, rho(jm), diff_m(jm), chie_m(jm), chii_m(jm)
      write(*,50) jm, rho(jm), exch_m(jm), egamma_m(jm), egamma_d(jm,1)
      write(*,50) jm, rho(jm), gamma_p_m(jm),anrate_m(jm),anrate2_m(jm)
      write(*,50) jm, rho(jm), anfreq_m(jm), anfreq2_m(jm)
      write(*,50) jm, rho(jm), gamma_p_exp(jm), vphi_m(jm)*corot,
     >            vpar_m(jm), vper_m(jm)
      write(*,50) jm, rho(jm), etaphi_m(jm), etapar_m(jm), etaper_m(jm)
c
      endif
c
      endif
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c... Simple chi=1 test (no neoclassical)
c
      if (imodel.eq.99) then
        diffnem=0.D0
        chietem=1.D0
        chiitim=1.D0
        chiegb_m(jm)=chietem/cgyrobohm_m(jm)
        chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
        diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
        use_xneo_m=-1
        chiineo_exp(jm)=0.D0
      endif
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c... Calculate neoclassical transport from model or experimental profiles
c
c... use_xneo_m switches for neoclassical transport
c       -1       use artificial chi near axis
c        0       use KAPISN w/ experimental profiles
c        1       use KAPISN w/ model profiles
c        2       use NCLASS w/ experimental profiles
c        3       use NCLASS w/ model profiles
c
c... Note: Adapted version KAPISN_ONE calculates chii at one grid point
c... Note that arrays are name.nc while single parameters are name.nco
c
c artificial chi at small r
c Note: use ineo=1
c
      if (use_xneo_m.eq.-1) then
        chiineo_m(jm)=chiineo_exp(jm)
        chieneo_m(jm)=chiineo_exp(jm)
        chiineogb_m(jm)=chiineo_m(jm)/cgyrobohm_m(jm)
        chieneogb_m(jm)=chieneo_m(jm)/cgyrobohm_m(jm)
c
      elseif (use_xneo_m.eq.0) then
        chiineo_m(jm)=chiineo_exp(jm)
        chiineogb_m(jm)=chiineogb_exp(jm)
c
      elseif (use_xneo_m.eq.1) then
c
        if (itest_exp.eq.7777) then
          tem=te_exp(jm)
          tim=ti_exp(jm)
          nem=ne_exp(jm)
          nim=ni_exp(jm)
        endif
c
        ng_nc=1                    !number of hyd. species
        aplasm_nc(1)=amassgas_exp  !array of atomic masses of hyd. species
        numzones_nc=nj_d           !number of zones
        btf_nc=bt_exp              !tf (tesla) fld at cntr of outer flux
        drshaf_nc=0.D0             !shaf shift of outer zone bdr. (cm)
c       rminor_nc=r_d(nj_d)*1.D2   !plasma minor radius (cm)
c       rmajor_nc=rmajor_exp*1.D2  !major radius (cntr of outer flux) (cm)
        rminor_nc=0.5D0*(rmin_exp(jm+1)+rmin_exp(jm))*1.D2
        rmajor_nc=rmaj_exp(jmaxm)*1.D2
c
        rhoel_nco=abs(nem*1.D13)        !electron density (cm**-3)
        rhob_nco(1)=abs(nim*1.D13)      !hyd. spec. den s (cm**-3)
        rhi_nco=abs(en_d(jm+1,2)*1.D-6) !z.c. array of av. impurity s density
        rhoi_nco=rhi_nco+rhob_nco(1)    !z.c. array of total ion density 
        te_nco=abs(tem*1.D3)            !z.c. array of Te (ev)
        ti_nco=abs(tim*1.D3)            !z.c. array of Ti (ev)
        zeff_nco=(zeff_exp(jm+1)+zeff_exp(jm))/2.D0  !z.c. array of plasma zeff
        q_nco=(q_exp(jm+1)+q_exp(jm))/2.D0           !z.c. array of safety factor
        aimp_nco=pimpa
        xzimp_nco=pimpz
c
c       write(*,'(i3,2x,0p1f4.2,1pe13.5,a10)') jm,
c    >       rho(jm), te_nco, ' model'
c        write(*,50) jm, rho(jm), tem, te_exp(jm)
c
c...    * means: default value given in input.mlt/mlt0in
c...    Note that the comments might be screwed up since some commentlines
c...    had to be deleted otherwise BASIS could not handle it (...) 
c...    See file neoclass.basf for the unscrewed comments
c
        call kapisn_one(
     >          nkimod_nc,         !*kapai model nr desired
     >          aimp_nc,           !*atomic mass of av. impurity
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
     >          rminor_nc,         !plasma minor radius (cm)
     >          rmajor_nc,         !major radius (cntr of outer flux) (cm)
     >          istringer_nc,      !Stringer correctioon
     >          xkapi_nco,         !o Neo-Class Ion thermal diff. (cm**2/sec)
     >          xnstari_nco,       !o nu-star-ions, 
     >          xnstare_nco,       !o nu-star-elecs, 
     >          ztau_nco,          !o ion collision time
     >          zrho_nco,          !o ion poloidal gyro-radius
     >          zfluxlim_nco)      !o flux lim flow max temp grad length
c
c       write(*,'(i3,2x,0p1f4.2,1pe13.5,a10)') jm, 
c     >       rho(jm), xkapi_nco, ' model'
       xkapi_nco=xkapi_nco*1.D-4*nim*1.D19            !1/(m*s)
     >     /elong_exp(jm)
c
       cgyrobohm_m(jm)=1.D-4*979.D3*(tem*1.D3)**0.5D0/
     >  (arho_exp*100.D0)*(102.D0*(tem*1.D3)**0.5D0/
     >  bt_exp/1.D4)**2.D0*(amassgas_exp)**0.5D0
c
c...   calculate conductivity

c
       chiineo_m(jm)=xkapi_nco/nim/1.D19/gradrhosq_exp(jm)    !m**2/s
       chiineogb_m(jm)=chiineo_m(jm)/cgyrobohm_m(jm)
c      write(*,'(i2,2x,0p1f4.2,1p4e13.5,a10)') jm, 
c    >       rho(jm), chiineo_exp(jm), chiineo_m(jm), 
c    >       q_exp(jm), q_d(jm),' model'
c
c      write(*,*) jm,' chieneogb_m = ',chiineogb_m(jm)
c
      elseif (use_xneo_m.eq.2) then  ! FORCEBAL neoclassical from exp.
        ineo=1  ! use chieneo in chietem and dineo in diffnem
        chieneo_m(jm)=chieneo_exp(jm)
        chiineo_m(jm)=chiineo_exp(jm)
        chieneogb_m(jm)=chieneogb_exp(jm)
        chiineogb_m(jm)=chiineogb_exp(jm)
        deneo_m(jm)=deneo_exp(jm)
        dineo_m(jm)=dineo_exp(jm)
c
      elseif (use_xneo_m.eq.3) then  ! call FORCEBAL from TRCOEF
        ineo=1  ! use chieneo in chietem and dineo in diffnem
        chieneo_m(jm)=0.D0
        chiineo_m(jm)=0.D0
        chieneogb_m(jm)=0.D0
        chiineogb_m(jm)=0.D0
        deneo_m(jm)=0.D0
        dineo_m(jm)=0.D0
      endif
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c  Add neoclassical transport to anomalous transport and compute totals
c
       fac_imp_flow=1.D0
       if(iexp_imp.eq.1) fac_imp_flow=1.D0+nz_exp(jm)/ni_exp(jm)
c
      powineo_m(jm)=fac_imp_flow*
     >  kevdsecpmw*tim*nim*1.D19/arho_exp*gradrhosq_exp(jm)*
     >       sfactor(jm)*(cneo*chiineo_m(jm)*zpmti)
c
      chiitim=chiitim+cneo*chiineo_m(jm)
c      if(i_proc.eq.0) then
c      write(*,50) jm, rho(jm), chiineo_m(jm), chiineo_exp(jm)
c      endif
c
      if (ineo.eq.1) then
         chietem=chietem+cneo*chieneo_m(jm)
      endif
      if(ineo.eq.0) chietem=chietem+
     >   cneo/sqrt(amassgas_exp)/42.D0*chiineo_m(jm)
      if(ineo.eq.-1) chietem=chietem+cneo*chiineo_m(jm)
      if(ineo.eq.-2) then
        if(q_exp(jm).lt.1.) then
          chietem=chietem+cneo*chiineo_m(jm)
        endif
        if(q_exp(jm).ge.1.) then
          chietem=chietem+
     >            cneo/sqrt(amassgas_exp)/42.D0*chiineo_m(jm)
        endif
      endif
      if(ineo.eq.-3) then
        if(q_exp(jm).lt.1.) chietem=chietem+
     >   cneo*chiineo_m(jm)
        if(q_exp(jm).ge.1.) chietem=chietem+
     >   cneo*chiineo_m(jm)
      endif
      if(ineo.eq.-4) then
        if(q_exp(jm).lt.1.) then
          chietem=chietem+5.0*cneo*chiineo_m(jm)
        endif
        if(q_exp(jm).ge.1.) then
          chietem=chietem+
     >            cneo/sqrt(amassgas_exp)/42.D0*chiineo_m(jm)
        endif
      endif
      if(ineo.eq.-12) then
        if(q_exp(jm).lt.1.0) chietem=chietem+
     >  cneo*chiineogb_m(jm)
        if(q_exp(jm).ge.1.0) chietem=chietem+
     >   10.D0*cneo*1.D0/sqrt(amassgas_exp)/42.D0*chiineogb_m(jm)
      endif
      if(ineo.eq.-10) chietem=chietem+
     > 10.D0*cneo/sqrt(amassgas_exp)/42.D0*chiineo_m(jm)
c
c add in neoclassical particle and momentum diffusivity
c Note: use chiineo for both if ineo.ne.1
c Use delta^3/2*chiineo in momentum diffusivity if ineophi=1
c Hinton, Wong, Phys Fluids 28, 3082 (1985).
c
      if (ineo.eq.1) then
        diffnem=diffnem+dineo_m(jm)
      else
        diffnem=diffnem+cneod/dsqrt(amassgas_exp)/
     >          42.D0*chiineo_m(jm)
        if (ineophi.eq.1) then
          etaphim=etaphim+
     >            cneophi*delta_exp(jm)**1.5D0*chiineo_m(jm)/6.6D0
        else
          etaphim=etaphim+cneophi*chiineo_m(jm)
        endif
      endif
      if(xparam_pt(11).gt.0) then
        etaphim = etaphim + xparam_pt(11)*DABS(chietem)
      endif
c     write(*,*) 'ineo = ',ineo
c     write(*,50) jm, rho(jm), delta_exp(jm), etaphim,
c    >            cneophi*delt_exp(jm)**1.5D0*chiineo_m(jm)/6.6D0,
c    >            chiitim
c    
c caution zpmni_q or zpmni in flow????
      zpmnix=zpmni
c
c add grid scale diffusion
c
c     write(*,*) 'adiffphi_dv = ',adiffphi_dv
      art_diff=adiff_dv*cgyrobohm_m(jm)*geofac(jm)
c     write(*,50) jm, rho(jm), dr(jm,2), zpmti,
c    >        chiitim,art_diff*DABS(dr(jm,2)*zpmti/arho_exp)
      diffnem=diffnem+art_diff*DABS(dr(jm,2)*zpmni/arho_exp)
      chiitim=chiitim+art_diff*DABS(dr(jm,2)*zpmti/arho_exp)
      chietem=chietem+art_diff*DABS(dr(jm,2)*zpmte/arho_exp)
      if (adiffphi_dv.eq.0) then
        etaphim=etaphim+
     >        art_diff*DABS(dr(jm,2)*gamma_p_m(jm)/arho_exp)
        etaperm=etaperm+
     >        art_diff*DABS(dr(jm,2)*gamma_p_m(jm)/arho_exp)
      else
        etaphim=etaphim+
     >       adiffphi_dv*DABS(dr(jm,2)*gamma_p_m(jm)/arho_exp)
      endif
c     write(*,50) jm, rho(jm), cgyrobohm_m(jm), etaphim
c
c add small amount of background diffusivity
c
      chiitim=chiitim+cartdiff
      chietem=chietem+cartdiff
      etaphim=etaphim+cartdiff
c
c total diffusivities
c
      chiegb_m(jm)=chietem/cgyrobohm_m(jm)
      chiigb_m(jm)=chiitim/cgyrobohm_m(jm)
      diffgb_m(jm)=diffnem/cgyrobohm_m(jm)
      etagb_phi_m(jm)=etaphim/cgyrobohm_m(jm)
c
      if(imodel.eq.8) zpmnix=zpmni_q
c
c compute fluxes
c
       nifluxm = nim*zpmni*diffnem/arho_exp
       tefluxm = tem*nem*zpmte*chietem/arho_exp
       tifluxm = tim*nim*zpmti*chiitim/arho_exp
       vphifluxm = -nim*etaphim*gamma_p_gf*csda_m(jm)/1000.D0
       vperfluxm = -nim*etaperm*gamma_e_gf*csda_m(jm)/1000.D0
c
       powem=kevdsecpmw*tem*nem*1.D19/arho_exp*gradrhosq_exp(jm)*
     >       sfactor(jm)*(chietem*zpmte+chietim*zpmti+chienem*zpmnix)
     > +xconv*1.5D0*tem*flow_exp(jm)
       powim=fac_imp_flow*
     >       kevdsecpmw*tim*nim*1.D19/arho_exp*gradrhosq_exp(jm)*
     >       sfactor(jm)*(chiitem*zpmte+chiitim*zpmti+chiinem*zpmnix)
     > +xconv*1.5D0*tim*flow_exp(jm)
       flowm=kevdsecpmw*1.*nem*1.D19/arho_exp*gradrhosq_exp(jm)*
     >       sfactor(jm)*(difftem*zpmte+difftim*zpmti+diffnem*zpmnix)
c
c       write(*,*) jm, sfactor(jm), powim, ' powim'
c       write(*,*) jm, nim, tim, ' nim,tim'
c       write(*,*) jm, zpmte, zpmti, ' zpmte, zpmti'
c       write(*,*) jm, zpte_m(jm), zpti_m(jm), ' zpte_m, zpti_m'
c       write(*,*) jm, fac_imp_flow, gradrhosq_exp(jm), ' fac'
c
crew       powe_m(jm)=powem
crem       powi_m(jm)=powim
crew       flow_m(jm)=flowm
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
 50    format(i2,2x,0p1f5.3,0p6e13.5)
 51    format(i2,2x,0p1f9.6,1p6e14.6)
 53    format('GLF23 diagnostic printout',/' grid pt = ',i3,', rho = ',
     >        0p1f5.3,', time-step = ',i5)
 55    format(' iky',3x,'ky',6x,'gamma_k',5x,'freq_k',5x,'phi-norm',4x,
     >        'chie-tot',4x,'chie-etg',4x,'chii-tot')
 60    format(i3,2x,0p1f6.4,0p8e12.4)
 65    format(/,' anrate = ',0p1e13.5,2x,'anfreq = ',0p1e13.5)
 70    format(' chie-glf = ',0p1e13.5,2x,'chii-glf = ',0p1e13.5,2x,
     >        'etaphi-glf = ',0p1e13.5)
 75    format(' chie_gf = ',0p1e13.5,2x,'chie-e-glf = ',0p1e13.5,2x,
     >        'chii-glf = ',0p1e13.5)
c
       return 
       end
