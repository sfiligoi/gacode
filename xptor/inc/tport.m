c@tport.m General Atomics
c 14-feb-11 removed chiegb_sol_m, chiigb_sol_m, diffgb_sol_m
c 11-feb-11 added volavebetat_m, volavebetae_m, volavebetai_m
c           and betat_m, betaf_m
c 25-aug-09 added curboot_exp, curbeam_exp, currf_exp, curohm_exp
c           sion_exp, srecom_exp, scx_exp, sbcx_exp, s_exp, sdot_exp
c           enn_exp, ennw_exp, ennv_exp, volsn_exp
c           sbeame_exp, sbeam_exp, torque_exp
c 21-feb-08 added er_exp
c 10-sep-07 added tefluxm_etg
c 20-feb-07 added nebar_exp
c 12-aug-06 removed nsm
c 25-aug-03 added nsm
c 10-jul-03 added dineo_exp, deneo_exp, vineo_exp, veneo_exp in common
c 26-jun-03 added veneo_m, vineo_m for neoclass. pinch velocities
c 21-nov-02 added pheat_exp for total heating power
c 12-may-02 added gtautot_exp for total taue from 0D file
c 02-may-02 added pfast_exp for fast ion pressure
c 25-feb-02 added sfarea_exp for flux surface area, volavene_exp, volaveti_m,te_m
c 21-feb-02 added zpti2_exp for 2nd derivative in Ti_exp
c 21-nov-01 added flow_recom_exp for source due to recombination
c 06-nov-01 added tauescale_exp(5) for various taue scalings
c 11-oct-01 added sbcx_exp(0:jmaxmt,2)
c 27-jul-01 added rlne_exp, rlni_exp, rlte_exp, rlti_exp
c 13-jul-01 added ineophi, delt_exp (use delta from kapisn for chi-phi)
c 05-may-01 added istep_smoo,smoo_tim, smoo_nim to smooth Tipro, dinpro in ptor.f
c 12-mar-01 added btscale to scale B_toroidal in rep_iter
c 29-jan-01 added snbscale to scale beam ion source, nscale for densities
c 25-jan-01 added powe_lh_exp for lower hybrid heating
c 19-jan-01 commented out arho_exp_s - fne_s variables since not used
c 18-jan-01 added winc_exp_s, winc_m_s, fw_s
c 18-jan-01 removed wem, wim, gradwem, gradwim
c           added gradtim, gradtem, fi_m, fz_m, vphim, fim
c 17-jan-01 added nrun_s, scanid, infile_s
c 08-jan-01 added wstr_exp, wstr_core_exp, wstr_ped_exp, 
c           wstr_m, wstr_core_m, wstr_ped_m
c 21-nov-00 added omexb_ncl_exp for omega_exb from NCLASS (rad/s)
c           also added egamma_ncl_exp
c 14-nov-00 added nm1_exp, nm2_exp
c 13-nov-00 added dineo_exp, deneo_exp, *_m, for neocl. ptcle transport
c 01-nov-00 added nz_t, nitot_exp, nzm, nz_m, zpmnz, zpnz_m, zpnz_exp,
c           zpmnitot, zpnitot_exp, zpnitot_m
c 03-oct-00 added powe_exp_t, powi_exp_t, powewdot_exp_t, powiwdot_exp_t
c           to store integrated power flows and write to NetCDF file
c 25-sept-00 added pbescale and pbiscale, prfscale
c 29-june-00 added vphi_ncl_exp and vpol_ncl_exp for NCLASS velocities
c 28-june-00 added chieneo_exp and chieneogb_exp
c 11-may-00 changed jmaxmt from 50 to 300
c           mask_r(0:50) -> mask_r(0:jmaxmt) and put before nstk_s(150)
c           other various reordering of common blocks due to increased
c           grid size.
c 11-feb-00 added vper_t, theta_exp, zptheta_exp, nimpm, zpmnimp
c           nifluxm,tifluxm,tefluxm,vphifluxm,vparfluxm,vperfluxm
c           wem,wim,gradwem,gradwim,gradnim,gradnem,gradvpolm
c           gradvphim,gradvparm,gradvperm,vperm, ntime_t
c 10-feb-00 added mscale, nfscale (for flow_exp)
c 18-jan-00 added egamma_vphi, egamma_vpol, egamma_vstar
c 5-jan-00 added variables to save diffusivities from GLF23
c 11-nov-99 increased kmaxmt from 1000 to 1500 and added 
c           vexb_t, vetor_t, vepol_t, vstar_t, egamma_t,
c           anrate_t, etaphi_t, and gammap_t
c 05-nov-99 added bteff_exp(0:jmaxmt) after csda_exp(0:jmaxmt)
c 29-sept-99 added NCLASS variables
c 11-sept-00 jek added te_exp_t, and ti_exp_t arrays
c 9-14-99 added real arrays:
c   neflux(0:jmaxmt), teflux(0:jmaxmt), niflux(0:jmaxmt)
c   tiflux(0:jmaxmt), vphiflux(0:jmaxmt)
c   vparflux(0:jmaxmt), vperflux(0:jmaxmy)
c 9-14-99 added integer array mask_r(0:50)
c   put nstk_s before mprint in common block
c---------------------------------------------------------------------
      integer jmaxmt, kmaxmt, stkmax, icalleimx,nspecies
      parameter (jmaxmt=300, kmaxmt=3800, stkmax=1, icalleimx=100)
      parameter (nspecies=3)
      complex*16 xi
c..301x2001 elements
      real*8 te_t(0:jmaxmt,0:kmaxmt), ti_t(0:jmaxmt,0:kmaxmt)
      real*8 te_exp_t(0:jmaxmt,0:kmaxmt), ti_exp_t(0:jmaxmt,0:kmaxmt)
      real*8 vphi_exp_t(0:jmaxmt,0:kmaxmt)
      real*8 ni_t(0:jmaxmt,0:kmaxmt), ne_t(0:jmaxmt,0:kmaxmt)
      real*8 nz_t(0:jmaxmt,0:kmaxmt), vpol_t(0:jmaxmt,0:kmaxmt)
      real*8 vexb_t(0:jmaxmt,0:kmaxmt), vphi_t(0:jmaxmt,0:kmaxmt)
      real*8 xte_t(0:jmaxmt,0:kmaxmt), xti_t(0:jmaxmt,0:kmaxmt)
      real*8 vetor_t(0:jmaxmt,0:kmaxmt)
      real*8 vepol_t(0:jmaxmt,0:kmaxmt), vstar_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_t(0:jmaxmt,0:kmaxmt),anratem_t(0:jmaxmt,0:kmaxmt)
      real*8 etapar_t(0:jmaxmt,0:kmaxmt),etaphi_t(0:jmaxmt,0:kmaxmt)
      real*8 gammap_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_vphi_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_vpol_t(0:jmaxmt,0:kmaxmt)
      real*8 egamma_vstar_t(0:jmaxmt,0:kmaxmt)
      real*8 chie_t(0:jmaxmt,0:kmaxmt), chii_t(0:jmaxmt,0:kmaxmt)
      real*8 pech_t(0:jmaxmt,0:kmaxmt)
      real*8 powe_t(0:jmaxmt,0:kmaxmt), powi_t(0:jmaxmt,0:kmaxmt)
c..3801 elements
      real*8 time_t(0:kmaxmt), pechtot_t(0:kmaxmt)
      real*8 powe_exp_t(0:kmaxmt), powi_exp_t(0:kmaxmt)
      real*8 powewdot_exp_t(0:kmaxmt), powiwdot_exp_t(0:kmaxmt)
c..3010 elements (was 510)
      real*8 egamma_d(0:jmaxmt,10)
c
c..1200 elements (was 200)
      real*8 powescan_m(4,1:jmaxmt), powiscan_m(4,1:jmaxmt)
      real*8 flowscan_m(4,1:jmaxmt),chiegbscan_m(4,1:jmaxmt)
      real*8 chiigbscan_m(4,1:jmaxmt),diffgbscan_m(4,1:jmaxmt)
      real*8 anratescan_m(4,1:jmaxmt),dnratescan_m(4,1:jmaxmt)
      real*8 dtnratescan_m(4,1:jmaxmt),anfreqscan_m(4,1:jmaxmt)
c..903 elements (was 151)
      real*8 flow_p_j(0:jmaxmt,3)
c..602 elements
      real*8 sion_exp(0:jmaxmt,2), srecom_exp(0:jmaxmt,2)
      real*8 scx_exp(0:jmaxmt,2), sbcx_exp(0:jmaxmt,2) 
      real*8 s_exp(0:jmaxmt,2), sdot_exp(0:jmaxmt,2)
      real*8 enn_exp(0:jmaxmt,2), ennw_exp(0:jmaxmt,2)
      real*8 ennv_exp(0:jmaxmt,2), volsn_exp(0:jmaxmt,2)
c
c..301 elements
      real*8 theta_exp(0:jmaxmt), zptheta_exp(0:jmaxmt)
c 301 elements (was 51)
      real*8 radfunc(0:jmaxmt), te_m(0:jmaxmt), ti_m(0:jmaxmt)
      real*8 fi_m(0:jmaxmt), fz_m(0:jmaxmt)
      real*8 ne_m(0:jmaxmt), ni_m(0:jmaxmt), ns_m(0:jmaxmt)
      real*8 nz_m(0:jmaxmt)
      real*8 interchange_DR_m(0:jmaxmt)
      real*8 neflux(0:jmaxmt), teflux(0:jmaxmt), niflux(0:jmaxmt)
      real*8 tiflux(0:jmaxmt), vphiflux(0:jmaxmt)
      real*8 vparflux(0:jmaxmt), vexbflux(0:jmaxmt)
      real*8 zpte_m(0:jmaxmt)
      real*8 zpti_m(0:jmaxmt), zpne_m(0:jmaxmt), zpni_m(0:jmaxmt)
      real*8 zpnz_m(0:jmaxmt), zpnitot_m(0:jmaxmt)
      real*8 zpti2_m(0:jmaxmt)
      real*8 te_exp(0:jmaxmt), ti_exp(0:jmaxmt), ne_exp(0:jmaxmt)
      real*8 tfast_exp(0:jmaxmt)
      real*8 te_exp_sav(0:jmaxmt), ti_exp_sav(0:jmaxmt)
      real*8 angrotp_exp_sav(0:jmaxmt)
      real*8 vphip_exp_sav(0:jmaxmt)
      real*8 ni_exp(0:jmaxmt), nz_exp(0:jmaxmt), nfast_exp(0:jmaxmt)
      real*8 ptot_exp(0:jmaxmt), pfast_exp(0:jmaxmt)
      real*8 torque_exp(0:jmaxmt)
      real*8 stress_tor_exp(0:jmaxmt),stress_par_exp(0:jmaxmt)
      real*8 nitot_exp(0:jmaxmt)
      real*8 nm1_exp(0:jmaxmt), nm2_exp(0:jmaxmt), nm3_exp(0:jmaxmt)
      real*8 vexb_exp(0:jmaxmt), vpol_exp(0:jmaxmt)
      real*8 vdia_exp(nspecies,0:jmaxmt), vneo_exp(nspecies,0:jmaxmt)
      real*8 mach_exp(nspecies,0:jmaxmt)
      real*8 vphi_exp(0:jmaxmt),vpar_exp(0:jmaxmt),vper_exp(0:jmaxmt)
      real*8 vphie_exp(0:jmaxmt),vpare_exp(0:jmaxmt)
      real*8 vphiz_exp(0:jmaxmt),vparz_exp(0:jmaxmt)
      real*8 pzmn_sol(0:jmaxmt), ne_p(0:jmaxmt)
      real*8 vphi_ncl_exp(0:jmaxmt), vpol_ncl_exp(0:jmaxmt)
      real*8 omexb_ncl_exp(0:jmaxmt)
      real*8 ni_p(0:jmaxmt), zeff_exp(0:jmaxmt), zpte_exp(0:jmaxmt)
      real*8 zpti_exp(0:jmaxmt), zpne_exp(0:jmaxmt), zpni_exp(0:jmaxmt)
      real*8 zpnz_exp(0:jmaxmt), zpnitot_exp(0:jmaxmt)
      real*8 zpne2_exp(0:jmaxmt)
      real*8 zpte2_exp(0:jmaxmt)
      real*8 zpti2_exp(0:jmaxmt)
      real*8 rlne_exp(0:jmaxmt), rlni_exp(0:jmaxmt)
      real*8 rlte_exp(0:jmaxmt), rlti_exp(0:jmaxmt)
      real*8 vol_exp(0:jmaxmt), volume_exp(0:jmaxmt) 
      real*8 powe_beam_exp(0:jmaxmt)
      real*8 powe_rf_exp(0:jmaxmt), powe_lh_exp(0:jmaxmt)
      real*8 powe_oh_exp(0:jmaxmt)
      real*8 powe_rad_exp(0:jmaxmt), powe_ion_exp(0:jmaxmt)
      real*8 powe_wdot_exp(0:jmaxmt), powe_fus_exp(0:jmaxmt)
      real*8 powi_beam_exp(0:jmaxmt), powi_rf_exp(0:jmaxmt)
      real*8 powi_ion_exp(0:jmaxmt), powi_cx_exp(0:jmaxmt)
      real*8 powi_wdot_exp(0:jmaxmt), powi_fus_exp(0:jmaxmt)
      real*8 pow_ei_exp(0:jmaxmt), pow_ei_cor_m(0:jmaxmt)
      real*8 pow_ei_m(0:jmaxmt), pow_ei_mexp(0:jmaxmt)
      real*8 pow_ei_glf(0:jmaxmt)
      real*8 exch_m(0:jmaxmt), exch_exp(0:jmaxmt)
      real*8 exch_glf(0:jmaxmt)
      real*8 powe_fus_m(0:jmaxmt), powi_fus_m(0:jmaxmt)
      real*8 powe_fus_cor_m(0:jmaxmt), powi_fus_cor_m(0:jmaxmt)
      real*8 pow_br_m(0:jmaxmt), pow_br_cor_m(0:jmaxmt)
      real*8 pow_cycl_m(0:jmaxmt), pow_cycl_cor_m(0:jmaxmt)
      real*8 flow_wall_exp(0:jmaxmt), flow_recom_exp(0:jmaxmt)
      real*8 flow_beam_exp(0:jmaxmt), flow_sdot_exp(0:jmaxmt)
      real*8 powineo_exp(0:jmaxmt)
      real*8 qbeame_exp(0:jmaxmt), qbeami_exp(0:jmaxmt)
      real*8 qrfe_exp(0:jmaxmt), qrfi_exp(0:jmaxmt)
      real*8 qohm_exp(0:jmaxmt), qrad_exp(0:jmaxmt)
      real*8 qione_exp(0:jmaxmt), qioni_exp(0:jmaxmt)
      real*8 sbeame_exp(0:jmaxmt), sbeam_exp(0:jmaxmt)
      real*8 powineo_m(0:jmaxmt), flow_exch_exp(0:jmaxmt)
      real*8 powe_exp(0:jmaxmt), powi_exp(0:jmaxmt),flow_exp(0:jmaxmt)
      real*8 powe_m(0:jmaxmt), powi_m(0:jmaxmt), powz_m(0:jmaxmt)
      real*8 flowe_m(0:jmaxmt),flowi_m(0:jmaxmt),flowz_m(0:jmaxmt)
      real*8 stress_tor_i_m(0:jmaxmt),stress_tor_z_m(0:jmaxmt)
      real*8 stress_par_i_m(0:jmaxmt),stress_par_z_m(0:jmaxmt)
      real*8 stress_par_cor_m(0:jmaxmt)
      real*8 powe_glf(0:jmaxmt),powi_glf(0:jmaxmt),powz_glf(0:jmaxmt)
      real*8 flowe_glf(0:jmaxmt),flowi_glf(0:jmaxmt),flowz_glf(0:jmaxmt)
      real*8 stress_tor_i_glf(0:jmaxmt),stress_tor_z_glf(0:jmaxmt)
      real*8 stress_par_i_glf(0:jmaxmt),stress_par_z_glf(0:jmaxmt)
      real*8 powe_neo(0:jmaxmt), powi_neo(0:jmaxmt), powz_neo(0:jmaxmt)
      real*8 flowe_neo(0:jmaxmt),flowi_neo(0:jmaxmt),flowz_neo(0:jmaxmt)
      real*8 stress_tor_i_neo(0:jmaxmt),stress_tor_z_neo(0:jmaxmt)
      real*8 stress_par_i_neo(0:jmaxmt),stress_par_z_neo(0:jmaxmt)
      real*8 powe_adhoc(0:jmaxmt), powi_adhoc(0:jmaxmt)
      real*8 powz_adhoc(0:jmaxmt)
      real*8 flowe_adhoc(0:jmaxmt),flowi_adhoc(0:jmaxmt)
      real*8 flowz_adhoc(0:jmaxmt)
      real*8 stress_tor_i_adhoc(0:jmaxmt),stress_tor_z_adhoc(0:jmaxmt)
      real*8 stress_par_i_adhoc(0:jmaxmt),stress_par_z_adhoc(0:jmaxmt)
      real*8 diff_exp(0:jmaxmt),chii_exp(0:jmaxmt)
      real*8 chie_exp(0:jmaxmt),eta_tor_exp(0:jmaxmt)
      real*8 egamma_g(0:jmaxmt)
      real*8 rhosda_m(0:jmaxmt), csda_m(0:jmaxmt)
      real*8 vexb_m(0:jmaxmt), vpol_m(0:jmaxmt)
      real*8 vphi_m(0:jmaxmt), vpar_m(0:jmaxmt), vper_m(0:jmaxmt)
      real*8 vphie_m(0:jmaxmt), vpare_m(0:jmaxmt)
      real*8 vphiz_m(0:jmaxmt), vparz_m(0:jmaxmt)
      real*8 vdia_m(nspecies,0:jmaxmt), vneo_m(nspecies,0:jmaxmt)
      real*8 vdia_new(nspecies,0:jmaxmt),vneo_new(nspecies,0:jmaxmt)
      real*8 mach_m(nspecies,0:jmaxmt)
      real*8 vmode_m(0:jmaxmt), vstar_m(0:jmaxmt), vstarp_m(0:jmaxmt)
      real*8 vetor_m(0:jmaxmt), vepol_m(0:jmaxmt), egamma_m(0:jmaxmt)
      real*8 egamma_vphi(0:jmaxmt), egamma_vpol(0:jmaxmt)
      real*8 egamma_vstar(0:jmaxmt)
      real*8 gamma_mode_m(0:jmaxmt), gamma_p_m(0:jmaxmt)
      real*8 cgyrobohm_m(0:jmaxmt), betae_m(0:jmaxmt)
      real*8 betai_m(0:jmaxmt), betat_m(0:jmaxmt), betaf_m(0:jmaxmt)
      real*8 alpha_m(0:jmaxmt), xnu_m(0:jmaxmt)
      real*8 nu_pol_m(0:jmaxmt),kpol_m(0:jmaxmt),nuei_m(0:jmaxmt)
      real*8 vnewstare_m(0:jmaxmt), vnewstari_m(0:jmaxmt)
      real*8 chiegb_m(0:jmaxmt), chiigb_m(0:jmaxmt)
      real*8 diffgb_m(0:jmaxmt), diffgb_im_m(0:jmaxmt)
      real*8 exchgb_m(0:jmaxmt), etagb_phi_m(0:jmaxmt)
      real*8 etagb_par_m(0:jmaxmt), etagb_per_m(0:jmaxmt)
      real*8 chiegb_glf(0:jmaxmt), chiigb_glf(0:jmaxmt)
      real*8 diffgb_glf(0:jmaxmt), etagb_phi_glf(0:jmaxmt)
      real*8 etagb_par_glf(0:jmaxmt), etagb_per_glf(0:jmaxmt)
      real*8 chie_e_gb_m(0:jmaxmt), chie_m(0:jmaxmt), chii_m(0:jmaxmt)
      real*8 diff_m(0:jmaxmt),chie_e_m(0:jmaxmt)
      real*8 etapar_m(0:jmaxmt), etaper_m(0:jmaxmt)
      real*8 etaexb_m(0:jmaxmt), etaphi_m(0:jmaxmt)
      real*8 rotstab_m(0:jmaxmt)
      real*8 diffmxlgb_m(0:jmaxmt), deltheta_m(0:jmaxmt)
      real*8 taue_m(0:jmaxmt), taui_m(0:jmaxmt), taut_m(0:jmaxmt)
      real*8 taupe_m(0:jmaxmt), taupi_m(0:jmaxmt), rhosda_exp(0:jmaxmt)
      real*8 csda_exp(0:jmaxmt), bteff_exp(0:jmaxmt)
      real*8 vstar_exp(0:jmaxmt)
      real*8 vexbc_exp(0:jmaxmt), vexbcc_exp(0:jmaxmt)
      real*8 vstarp_exp(0:jmaxmt), egamma_exp(0:jmaxmt)
      real*8 egamma_ncl_exp(0:jmaxmt)
      real*8 megamma_exp(0:jmaxmt), gamma_p_exp(0:jmaxmt)
      real*8 gamma_mode_exp(0:jmaxmt), anrate_exp(0:jmaxmt)
      real*8 anfreq_exp(0:jmaxmt), cgyrobohm_exp(0:jmaxmt)
      real*8 betae_exp(0:jmaxmt), betai_exp(0:jmaxmt)
      real*8 alpha_exp(0:jmaxmt), xnu_exp(0:jmaxmt)
      real*8 vnewstare_exp(0:jmaxmt), vnewstari_exp(0:jmaxmt)
      real*8 chiegb_exp(0:jmaxmt), chiigb_exp(0:jmaxmt)
      real*8 diffgb_exp(0:jmaxmt), chiineogb_exp(0:jmaxmt)
      real*8 chiineo_exp(0:jmaxmt), chiineogb_m(0:jmaxmt)
      real*8 chiineo_m(0:jmaxmt), chiicond_exp(0:jmaxmt)
      real*8 chieneo_exp(0:jmaxmt), chieneogb_exp(0:jmaxmt)
      real*8 chieneogb_m(0:jmaxmt), chieneo_m(0:jmaxmt)
      real*8 chiecond_exp(0:jmaxmt)
      real*8 dineo_exp(0:jmaxmt), deneo_exp(0:jmaxmt)
      real*8 dineo_m(0:jmaxmt), deneo_m(0:jmaxmt)
      real*8 vineo_m(0:jmaxmt), veneo_m(0:jmaxmt)
      real*8 vineo_exp(0:jmaxmt), veneo_exp(0:jmaxmt)
      real*8 etaphineo_m(0:jmaxmt), etaperneo_m(0:jmaxmt)
      real*8 taue_exp(0:jmaxmt), taui_exp(0:jmaxmt), taut_exp(0:jmaxmt)
      real*8 tautot_exp(0:jmaxmt)
      real*8 taupe_exp(0:jmaxmt), taupi_exp(0:jmaxmt), rho(0:jmaxmt)
      real*8 rmin_exp(0:jmaxmt), rmaj_exp(0:jmaxmt), sfactor(0:jmaxmt)
      real*8 geofac(0:jmaxmt), drhodr(0:jmaxmt), drhodrrrho(0:jmaxmt)
      real*8 georotrate(0:jmaxmt), geoalpha(0:jmaxmt)
      real*8 h_exp(0:jmaxmt), f_exp(0:jmaxmt), g_exp(0:jmaxmt)
      real*8 bp_exp(0:jmaxmt)
      real*8 c_par(0:jmaxmt), c_per(0:jmaxmt), c_tor(0:jmaxmt)
      real*8 a_pol(0:jmaxmt), a_tor(0:jmaxmt), Bp0(0:jmaxmt)
      real*8 xb2_exp(0:jmaxmt), xbm2_exp(0:jmaxmt), xngrth_exp(0:jmaxmt)
     & ,     xgrbm2_exp(0:jmaxmt), fm1_exp(0:jmaxmt), fm2_exp(0:jmaxmt)
     & ,     fm3_exp(0:jmaxmt), fhat_exp(0:jmaxmt), xr2_exp(0:jmaxmt)
      real*8 gradrhosq_exp(0:jmaxmt), gradrho_exp(0:jmaxmt)
      real*8 sfarea_exp(0:jmaxmt), cxarea_exp(0:jmaxmt)
      real*8 psir_exp(0:jmaxmt)
      real*8 elong_exp(0:jmaxmt), delta_exp(0:jmaxmt)
      real*8 delt_exp(0:jmaxmt), q_exp(0:jmaxmt)
      real*8 curden_exp(0:jmaxmt), curboot_exp(0:jmaxmt) 
      real*8 curbeam_exp(0:jmaxmt), currf_exp(0:jmaxmt)
      real*8 curohm_exp(0:jmaxmt)
      real*8 shat_exp(0:jmaxmt), er_exp(0:jmaxmt)
      real*8 angrot_exp(0:jmaxmt), angrotp_exp(0:jmaxmt)
      real*8 psi_exp(0:jmaxmt), betat_exp(0:jmaxmt)
      real*8 b_unit_loc_s(0:jmaxmt), b_norm_loc_s(0:jmaxmt)
      real*8 rhosda_loc_s(0:jmaxmt), csda_loc_s(0:jmaxmt)
      real*8 aspectratio_loc_s(0:jmaxmt), kappa_loc_s(0:jmaxmt)
      real*8 s_kappa_loc_s(0:jmaxmt), delta_loc_s(0:jmaxmt)
      real*8 s_delta_loc_s(0:jmaxmt), shift_loc_s(0:jmaxmt)
      real*8 beta_loc_s(0:jmaxmt), q_loc_s(0:jmaxmt)
      real*8 shat_loc_s(0:jmaxmt), dlntidr_loc_s(0:jmaxmt)
      real*8 dlntedr_loc_s(0:jmaxmt), dlnnidr_loc_s(0:jmaxmt)
      real*8 dlnnedr_loc_s(0:jmaxmt), tiote_loc_s(0:jmaxmt)
      real*8 nione_loc_s(0:jmaxmt), xnu_loc_s(0:jmaxmt)
      real*8 shat_mhd_loc_s(0:jmaxmt), alpha_mhd_loc_s(0:jmaxmt)
      real*8 alpha_mhd1(0:jmaxmt), alpha_mhd2(0:jmaxmt)
      real*8 beta_loc_0_s(0:jmaxmt), volume_loc_s(0:jmaxmt)
      real*8 delta_dor_s(0:jmaxmt), d_prime_dor_s(0:jmaxmt)
      real*8 k_prime_dor_s(0:jmaxmt)
c
c 150 elements
c      real*8 arho_exp_s(stkmax), rmajor_exp_s(stkmax)
c      real*8 bt_exp_s(stkmax), tocur_d_s(stkmax)
c      real*8 elonga_exp_s(stkmax), deltaa_exp_s(stkmax)
c      real*8 amassgas_exp_s(stkmax), volaveniti_m_s(stkmax)
c      real*8 gtaut_m_s(stkmax), gtaui_m_s(stkmax)
c      real*8 n_t_tau_m_s(stkmax), n0_t0_tau_m_s(stkmax)
c      real*8 volaveniti2_m_s(stkmax), wstre_m_s(stkmax)
c      real*8 wstri_m_s(stkmax), te20_m_s(stkmax)
c      real*8 ti20_m_s(stkmax), tauiter_exp_s(stkmax)
c      real*8 p_glob_exp_s(stkmax), w_glob_exp_s(stkmax)
c      real*8 volaveniti_exp_s(stkmax), gtaut_exp_s(stkmax)
c      real*8 gtaui_exp_s(stkmax), n_t_tau_exp_s(stkmax)
c      real*8 n0_t0_tau_exp_s(stkmax), volaveniti2_exp_s(stkmax)
c      real*8 wstre_exp_s(stkmax), wstri_exp_s(stkmax)
c      real*8 wstre_ped_exp_s(stkmax), wstri_ped_exp_s(stkmax)
c      real*8 rhosda20_exp_s(stkmax), te20_exp_s(stkmax)
c      real*8 ti20_exp_s(stkmax), te_mix_s(stkmax)
c      real*8 ti_mix_s(stkmax), chisqte_s(stkmax)
c      real*8 chisqti_s(stkmax), chisqne_s(stkmax)
c      real*8 sigmate_s(stkmax), sigmati_s(stkmax)
c      real*8 sigmane_s(stkmax), fte_s(stkmax)
c      real*8 fti_s(stkmax), fne_s(stkmax)
c
c 100 elements
      real*8 te_cen_ei(1:icalleimx), te_mid_ei(1:icalleimx)
c
c 50 elements
      real*8 vscan_m(1:50)
c 9 elements
      real*8 dflowdzp(3,3),rmatrix(3,3),one(3,3)
c 5 elements
      real*8 tauescale_exp(5)
c 3 elements
      real*8 flow_p(3), zm(3), zpl(3), zph(3)
      real*8 zpsol(3),flowl(3),flowexp(3),fvector(3)
      real*8 vdia(nspecies),vneo(nspecies)
c
      real*8 cv,csdam
      real*8 pi_m,kevdsecpmw,tem,tim,nem,nim,nimpm,nzm
      real*8 zmte,zmti,zmne,zmni,zpmte,zpmnimp
      real*8 zpmti,zpmne,zpmni,zpmnz,zpmnitot
      real*8 pzmn(3),pzmn1,pzmn2,powem,powim
      real*8 flowm,chietem,chiitim,diffnem
      real*8 etaparm,etaperm,etaexbm,etaphim
      real*8 detamatrix,relxei,relxbr
      real*8 relxfus,tau_he,f_defect,f_def_beryllium,relx
      real*8 eps_tetest,ascale,bscale,pscale,mscale,nfscale
      real*8 te_mix,ti_mix
      real*8 chisqte,chisqti,chisqne
      real*8 pbescale, pbiscale, snbscale
      real*8 prfscale, prfescale, prfiscale
      real*8 nscale, btscale
      real*8 sigmate,sigmati,sigmane,fte,fti,fne
      real*8 neline_exp,niline_exp
      real*8 zpmte1, zpmti1,zpmne1,zpmni1,volaveniti_m
      real*8 volaveti_m, volavete_m
      real*8 gtaut_m,gtaui_m,n_t_tau_m
      real*8 n0_t0_tau_m,volaveniti2_m,wstre_m
      real*8 wstri_m,te20_m,ti20_m, volavene_exp
      real*8 volaveniti_exp,gtaut_exp,gtaui_exp
      real*8 gtauth_exp,gtautot_exp,pheat_exp
      real*8 n_t_tau_exp,n0_t0_tau_exp,volaveniti2_exp
      real*8 wstre_exp,wstri_exp,wstre_ped_exp,wstri_ped_exp
      real*8 wstr_exp, wstr_inc_exp, wstr_core_exp, wstr_ped_exp
      real*8 wstr_m, wstr_inc_m, wstr_core_m
      real*8 volavebetat_m, volavebetae_m, volavebetai_m
      real*8 volavebetaf_m, p_glob_exp,w_glob_exp
      real*8 rhosda20_exp,te20_exp,ti20_exp,nebar_exp
      real*8 rho_g,btvac_toq,btplas_toq
      real*8 torfO2pi_toq,a_bnd_toq,a_unit_exp,bt_mag_center
      real*8 alpha_e_dk,alpha_star_dk
      real*8 xp_time
      real*8 neflux_glf,niflux_glf,nzflux_glf
      real*8 teflux_glf,tiflux_glf,tzflux_glf
      real*8 vphiflux_glf,vphizflux_glf
      real*8 vparflux_glf,vparzflux_glf
      real*8 neflux_neo,niflux_neo,nzflux_neo
      real*8 teflux_neo,tiflux_neo,tzflux_neo
      real*8 vphiflux_neo,vphizflux_neo
      real*8 vparflux_neo,vparzflux_neo
      real*8 neflux_adhoc,niflux_adhoc,nzflux_adhoc
      real*8 teflux_adhoc,tiflux_adhoc,tzflux_adhoc
      real*8 vphiflux_adhoc,vphizflux_adhoc
      real*8 vparflux_adhoc,vparzflux_adhoc
      real*8 nefluxm,nifluxm,nzfluxm
      real*8 tifluxm,tefluxm,vparfluxm,vexbfluxm
      real*8 tefluxm_etg,vphifluxm,vperfluxm
      real*8 vphizfluxm,vparzfluxm,tzfluxm
      real*8 gradnim,gradnem,gradnzm,gradtim,gradtem
      real*8 fim,fzm,gradfim,gradfzm
      real*8 gradvexbm,gradvpolm,vexbm,vpolm
      real*8 smoo_tim, smoo_nim, smoo_vphim
      real*8 sum_ne_exp,sum_ne_m,ne_gain_pt,wall_mult
      real*8 sum_te_exp,sum_te_m,te_gain_pt,te_bc_mult
      real*8 sum_ti_exp,sum_ti_m,ti_gain_pt,ti_bc_mult
c
      integer mask_r(0:jmaxmt), nstk_s(150), mprint(10)
      integer iskip(stkmax), istk_smooth(stkmax)
      integer jstmix_s(stkmax)
      integer iei,iquit,iquitmax,ibifurcate,imodcall,jm,iloop_p_max
      integer iterate,iscale,iexp_exch,iexp_q,iexp_neo
      integer iexp_imp,ineo,ineophi
      integer jstmix,jstinv,ismooth,istep_smoo
      integer i_delay,i_gridsmooth,iv_dk,itest_dk, kmaxm, ntime_t
      integer dilution_model
c
      character*50 cudir
      character*(6) phase
      character*(10) tok
      character*(40) shot
c
      common /tport/ xi, te_t, ti_t, te_exp_t, ti_exp_t, vphi_exp_t
     & , ne_t, ni_t, nz_t, powe_t, powi_t, vphi_t, vexb_t, vpol_t
     & , xte_t, xti_t, vetor_t, vepol_t, vstar_t, egamma_t 
     & , egamma_vphi_t, egamma_vpol_t, egamma_vstar_t, anratem_t
     & , chie_t, chii_t, etapar_t, etaphi_t
     & , gammap_t, pech_t, time_t
     & , pechtot_t, powe_exp_t, powi_exp_t 
     & , powewdot_exp_t, powiwdot_exp_t, egamma_d
     & , powescan_m, powiscan_m, flowscan_m, chiegbscan_m
     & , chiigbscan_m,diffgbscan_m, anratescan_m,dnratescan_m
     & , dtnratescan_m,anfreqscan_m, flow_p_j
     & , sion_exp, srecom_exp, scx_exp, sbcx_exp, s_exp, sdot_exp
     & , enn_exp, ennw_exp, ennv_exp, volsn_exp
     & , radfunc, te_m, ti_m, ne_m, ni_m, ns_m, nz_m
     & , fi_m, fz_m
     & , interchange_DR_m
     & , neflux, teflux, niflux
     & , tiflux, vphiflux, vparflux, vexbflux
     & , zpte_m, zpti_m
     & , zpne_m, zpni_m, zpnz_m, zpnitot_m, zpti2_m
     & , te_exp, ti_exp, tfast_exp
     & , ne_exp, ni_exp, nz_exp, nfast_exp, nitot_exp
     & , nm1_exp, nm2_exp, nm3_exp
     & , te_exp_sav, ti_exp_sav, angrotp_exp_sav
     & , ptot_exp, pfast_exp, torque_exp
     & , stress_tor_exp, stress_par_exp, stress_par_cor_m
     & , flowe_glf, flowi_glf, flowz_glf
     & , powe_glf, powi_glf, powz_glf
     & , stress_tor_i_glf, stress_tor_z_glf
     & , stress_par_i_glf, stress_par_z_glf
     & , flowe_neo, flowi_neo, flowz_neo
     & , powe_neo, powi_neo, powz_neo
     & , stress_tor_i_neo, stress_tor_z_neo
     & , stress_par_i_neo, stress_par_z_neo
     & , flowe_adhoc, flowi_adhoc, flowz_adhoc
     & , powe_adhoc, powi_adhoc, powz_adhoc
     & , stress_tor_i_adhoc, stress_tor_z_adhoc
     & , stress_par_i_adhoc, stress_par_z_adhoc
     & , vpar_exp, vexb_exp, vper_exp
     & , vphie_exp, vpare_exp, vphiz_exp, vparz_exp
     & , vpol_exp, vdia_exp, vneo_exp, mach_exp
     & , vphi_exp, vphi_ncl_exp, vpol_ncl_exp
     & , omexb_ncl_exp, pzmn_sol, ne_p
     & , ni_p, zeff_exp, zpte_exp, zpti_exp
     & , zpne_exp, zpni_exp, zpnitot_exp, zpti2_exp
     & , zpte2_exp, zpne2_exp
     & , rlne_exp, rlni_exp, rlte_exp, rlti_exp
     & , zpnz_exp, vol_exp, volume_exp, powe_beam_exp, powe_rf_exp
     & , powe_lh_exp, powe_oh_exp
     & , powe_rad_exp, powe_ion_exp, powe_wdot_exp, powe_fus_exp
     & , powi_beam_exp, powi_rf_exp, powi_ion_exp, powi_cx_exp
     & , powi_wdot_exp, powi_fus_exp, pow_ei_exp, pow_ei_cor_m
     & , pow_ei_m, pow_ei_glf, pow_ei_mexp
     & , exch_m, exch_exp, exch_glf
     & , powe_fus_m, powi_fus_m, powe_fus_cor_m, powi_fus_cor_m
     & , pow_br_m, pow_br_cor_m, pow_cycl_m, pow_cycl_cor_m
     & , flow_wall_exp, flow_recom_exp, flow_beam_exp
     & , flow_sdot_exp, powineo_exp
     & , powineo_m, flow_exch_exp, powe_exp, powi_exp, flow_exp
     & , qbeame_exp, qbeami_exp, qrfe_exp, qrfi_exp, qohm_exp
     & , qrad_exp, qione_exp, qioni_exp, sbeame_exp, sbeam_exp
     & , powe_m, powi_m, powz_m
     & , flowe_m, flowi_m, flowz_m
     & , stress_tor_i_m, stress_tor_z_m
     & , stress_par_i_m, stress_par_z_m
     & , diff_exp, chie_exp, chii_exp, eta_tor_exp, egamma_g
     & , rhosda_m, csda_m, vexb_m, vpar_m, vphi_m, vper_m
     & , vphie_m, vpare_m, vphiz_m, vparz_m
     & , vdia_m, vneo_m, vdia_new, vneo_new
     & , vpol_m, nu_pol_m, kpol_m, nuei_m
     & , vmode_m, vstar_m, vstarp_m, vetor_m, vepol_m, egamma_m
     & , mach_m, egamma_vphi, egamma_vpol, egamma_vstar
     & , gamma_mode_m, gamma_p_m, cgyrobohm_m
     & , betae_m, betai_m, betat_m, betaf_m
     & , alpha_m, xnu_m, vnewstare_m, vnewstari_m, chiegb_m, chiigb_m
     & , diffgb_m, diffgb_im_m, exchgb_m, etagb_phi_m
     & , etagb_par_m, etagb_per_m
     & , chiegb_glf, chiigb_glf, diffgb_glf, etagb_phi_glf
     & , etagb_par_glf, etagb_per_glf
     & , chie_e_gb_m, chie_m, chii_m, diff_m, chie_e_m
     & , etapar_m, etaper_m, etaexb_m, etaphi_m
     & , rotstab_m, diffmxlgb_m, deltheta_m
     & , taue_m, taui_m, taut_m, taupe_m, taupi_m, rhosda_exp
     & , csda_exp, bteff_exp, vstar_exp, vexbc_exp
     & , vexbcc_exp, vstarp_exp, egamma_exp, egamma_ncl_exp
     & , megamma_exp, gamma_p_exp
     & , gamma_mode_exp, anrate_exp, anfreq_exp, cgyrobohm_exp
     & , betae_exp, betai_exp, alpha_exp, xnu_exp
     & , vnewstare_exp, vnewstari_exp, chiegb_exp, chiigb_exp
     & , diffgb_exp, chiineogb_exp, chiineo_exp
     & , chieneo_exp, chieneogb_exp, veneo_exp, vineo_exp
     & , deneo_exp, dineo_exp
     & , chiineogb_m, chiineo_m, chieneogb_m, chieneo_m
     & , deneo_m, dineo_m, veneo_m, vineo_m
     & , etaphineo_m, etaperneo_m
     & , chiicond_exp, chiecond_exp, taue_exp
     & , taui_exp, taut_exp, tautot_exp, taupe_exp, taupi_exp, rho
     & , rmin_exp, rmaj_exp, sfactor, geofac, drhodr, drhodrrrho
     & , georotrate, geoalpha
     & , h_exp, f_exp, g_exp,  psir_exp
     & , bp_exp, c_par, c_per, c_tor, a_pol, a_tor, Bp0
     & , xb2_exp, xr2_exp, xbm2_exp
     & , xngrth_exp, xgrbm2_exp,fm1_exp, fm2_exp, fm3_exp
     & , fhat_exp, gradrhosq_exp, gradrho_exp
     & , sfarea_exp, cxarea_exp
     & , elong_exp, delta_exp, delt_exp, q_exp, shat_exp
     & , curden_exp, curboot_exp, curbeam_exp
     & , currf_exp, curohm_exp
     & , er_exp, angrot_exp, angrotp_exp, psi_exp, betat_exp
     & , theta_exp, zptheta_exp
     & , b_unit_loc_s, b_norm_loc_s, rhosda_loc_s, csda_loc_s
     & , aspectratio_loc_s, kappa_loc_s, s_kappa_loc_s, delta_loc_s
     & , s_delta_loc_s, shift_loc_s, beta_loc_s, q_loc_s
     & , shat_loc_s, dlntidr_loc_s, dlntedr_loc_s, dlnnidr_loc_s
     & , dlnnedr_loc_s, tiote_loc_s, nione_loc_s, xnu_loc_s
     & , shat_mhd_loc_s, alpha_mhd_loc_s, alpha_mhd1, alpha_mhd2
     & , beta_loc_0_s, volume_loc_s, delta_dor_s, d_prime_dor_s
     & , k_prime_dor_s
c     & , arho_exp_s, rmajor_exp_s, bt_exp_s, tocur_d_s
c     & , elonga_exp_s, deltaa_exp_s, amassgas_exp_s, volaveniti_m_s
c     & , gtaut_m_s, gtaui_m_s, n_t_tau_m_s, n0_t0_tau_m_s
c     & , volaveniti2_m_s, wstre_m_s, wstri_m_s, te20_m_s
c     & , ti20_m_s, tauiter_exp_s,  p_glob_exp_s, w_glob_exp_s
c     & , volaveniti_exp_s, gtaut_exp_s, gtaui_exp_s, n_t_tau_exp_s
c     & , n0_t0_tau_exp_s, volaveniti2_exp_s, wstre_exp_s, wstri_exp_s
c     & , wstre_ped_exp_s, wstri_ped_exp_s
c     & , rhosda20_exp_s, te20_exp_s
c     & , ti20_exp_s, te_mix_s, ti_mix_s, chisqte_s
c     & , chisqti_s, chisqne_s, sigmate_s, sigmati_s
c     & , sigmane_s, fte_s, fti_s, fne_s
       common /tport1/ te_cen_ei, te_mid_ei, vscan_m
     & , dflowdzp, rmatrix, one, tauescale_exp
     & , flow_p, zm, zpl, zph, zpsol,flowl,flowexp,fvector
     & , vdia,vneo
     & , cv,csdam
     & , pi_m,kevdsecpmw,tem,tim,nem,nim,nimpm,nzm
     & , zmte,zmti,zmne,zmni
     & , zpmte,zpmti,zpmne,zpmni,zpmnz,zpmnitot,zpmnimp,pzmn
     & , pzmn1,pzmn2,powem,powim,flowm,chietem,chiitim,diffnem
     & , etaparm,etaperm,etaexbm,etaphim
     & , detamatrix,relxei,relxbr
     & , relxfus,tau_he,f_defect,f_def_beryllium,relx,eps_tetest
     & , ascale,bscale,pscale,mscale,nfscale,te_mix,ti_mix
     & , pbescale,pbiscale,snbscale,prfscale,prfescale,prfiscale
     & , nscale,btscale,chisqte,chisqti
     & , chisqne,sigmate,sigmati,sigmane,fte,fti,fne
     & , neline_exp,niline_exp,zpmte1
     & , zpmti1,zpmne1,zpmni1,volaveniti_m,volaveti_m,volavete_m
     & , gtaut_m,gtaui_m,n_t_tau_m,n0_t0_tau_m,volaveniti2_m,wstre_m
     & , wstri_m,te20_m,ti20_m,volavene_exp,volaveniti_exp
     & , gtaut_exp,gtaui_exp,gtauth_exp,gtautot_exp,pheat_exp
     & , n_t_tau_exp,n0_t0_tau_exp,volaveniti2_exp,wstre_exp,wstri_exp
     & , wstre_ped_exp,wstri_ped_exp,p_glob_exp,w_glob_exp
     & , rhosda20_exp,te20_exp,ti20_exp,nebar_exp
     & , rho_g,btvac_toq,btplas_toq
     & , torfO2pi_toq,a_bnd_toq,a_unit_exp,bt_mag_center,alpha_e_dk
     & , alpha_star_dk,xp_time
     & , neflux_glf,niflux_glf,nzflux_glf
     & , teflux_glf,tiflux_glf,tzflux_glf
     & , vphiflux_glf,vphizflux_glf
     & , vparflux_glf,vparzflux_glf
     & , neflux_neo,niflux_neo,nzflux_neo
     & , teflux_neo,tiflux_neo,tzflux_neo
     & , vphiflux_neo,vphizflux_neo
     & , vparflux_neo,vparzflux_neo
     & , neflux_adhoc,niflux_adhoc,nzflux_adhoc
     & , teflux_adhoc,tiflux_adhoc,tzflux_adhoc
     & , vphiflux_adhoc,vphizflux_adhoc
     & , vparflux_adhoc,vparzflux_adhoc
     & , nefluxm,nifluxm,nzfluxm
     & , tifluxm,tefluxm,vparfluxm,vexbfluxm
     & , tefluxm_etg,vphifluxm,vperfluxm
     & , vphizfluxm,vparzfluxm, tzfluxm
     & , gradnim,gradnem,gradnzm,gradtim,gradtem
     & , fim,fzm,gradfim,gradfzm
     & , gradvexbm,gradvpolm,vexbm,vpolm
     & , smoo_tim,smoo_nim,smoo_vphim
     & , sum_ne_exp,sum_ne_m,ne_gain_pt,wall_mult
     & , sum_te_exp,sum_te_m,te_gain_pt,te_bc_mult
     & , sum_ti_exp,sum_ti_m,ti_gain_pt,ti_bc_mult
     & , wstr_m, wstr_inc_m, wstr_core_m, wstr_exp
     & , wstr_inc_exp, wstr_core_exp, wstr_ped_exp
     & , volavebetat_m, volavebetae_m, volavebetai_m, volavebetaf_m
     & , mask_r, nstk_s, mprint, iskip, istk_smooth, jstmix_s
     & , iei,iquit,iquitmax,ibifurcate,imodcall,jm,iloop_p_max
     & , iterate,iscale,iexp_exch,iexp_q,iexp_neo
     & , iexp_imp,ineo,ineophi,dilution_model
     & , jstmix,jstinv,ismooth,istep_smoo
     & , i_delay,i_gridsmooth,iv_dk,itest_dk,kmaxm,ntime_t
     & , cudir, tok, shot, phase


