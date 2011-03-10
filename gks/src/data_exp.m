c@data_exp.m 30-may-02 G. Staebler, General Atomics
c common block for _exp variables for GKS stand alone code
c 25-sep-11 added torque_exp, sfarea_exp, cxarea_exp, psir_exp
c           f_exp, g_exp, volume_exp, curboot_exp, curbeam_exp, currf_exp
c           curohm_exp, qbeame_exp, qbeami_exp, qohm_exp, qrad_exp, 
c           qione_exp, qioni_exp, sbeame_exp, sbeam_exp, sion_exp
c           srecom_exp, s_exp, sdot_exp, enn_exp, ennw_exp, ennv_exp
c           volsn_exp, er_exp, prfescale, prfiscale, nebar_exp
c           removed btscale
c 12-may-02 added gtautot_exp for total taue from 0D file
c 02-may-02 added pfast_exp for fast ion pressure
c 25-feb-02 added vsurf_exp for flux surface area, volavene_exp, volaveti_m,te_m
c 21-feb-02 added zpti2_exp for 2nd derivative in Ti_exp
c 21-nov-01 added flow_recom_exp for source due to recombination
c 06-nov-01 added tauescale_exp(5) for various taue scalings
c 11-oct-01 added sbcx_exp(0:jmaxmt,2)
c 27-jul-01 added rlne_exp, rlni_exp, rlte_exp, rlti_exp
c 20-jul-01 added idiffphi_dv to allow adiff to not be used with chi-phi
c 13-jul-01 added ineophi, delt_exp (use delta from kapisn for chi-phi)
c 05-may-01 added istep_smoo,smoo_tim, smoo_nim to smooth Tipro, dinpro in ptor.f
c 12-mar-01 added btscale to scale B_toroidal in rep_iter
c 29-jan-01 added snbscale to scale beam ion source, nscale for densities
c 25-jan-01 added powe_lh_exp for lower hybrid heating
c           added gradtim, gradtem, fi_m, fz_m, vphim, fim
c 17-jan-01 added nrun_s, scanid, infile_s
c 08-jan-01 added wstr_exp, wstr_core_exp, wstr_ped_exp
c 21-nov-00 added omexb_ncl_exp for omega_exb from NCLASS (rad/s)
c           also added egamma_ncl_exp
c 14-nov-00 added nm1_exp, nm2_exp
c 13-nov-00 added dineo_exp, deneo_exp, for neocl. ptcle transport
c 01-nov-00 added nz_t, nitot_exp, nzm, zpmnz,  zpnz_exp,
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
c           wem,wim,gradwem,gradwim,gradnim,gradnem
c           gradvphim,gradvparm,gradvperm,vperm, ntime_t
c 10-feb-00 added mscale, nfscale (for flow_exp)
c 18-jan-00 added egamma_vphi, egamma_vpol, egamma_vstar
c 5-jan-00 added variables to save diffusivities from GLF23
c 05-nov-99 added bteff_exp(0:jmaxmt) after csda_exp(0:jmaxmt)
c   put nstk_s before mprint in common block
c---------------------------------------------------------------------
      integer jmaxmt, kmaxmt, stkmax, icalleimx
      parameter (jmaxmt=200, kmaxmt=10, stkmax=1, icalleimx=1)
c      complex*16 xi
c
c..602 elements
      real*8 sion_exp(0:jmaxmt,2), srecom_exp(0:jmaxmt,2)
      real*8 scx_exp(0:jmaxmt,2), sbcx_exp(0:jmaxmt,2)
      real*8 s_exp(0:jmaxmt,2), sdot_exp(0:jmaxmt,2)
      real*8 enn_exp(0:jmaxmt,2), ennw_exp(0:jmaxmt,2)
      real*8 ennv_exp(0:jmaxmt,2), volsn_exp(0:jmaxmt,2)
c
c 201 elements (was 51)
      real*8 te_exp(0:jmaxmt), ti_exp(0:jmaxmt), ne_exp(0:jmaxmt)
      real*8 ni_exp(0:jmaxmt), nz_exp(0:jmaxmt), nfst_exp(0:jmaxmt)
      real*8 nitot_exp(0:jmaxmt)
      real*8 nm1_exp(0:jmaxmt), nm2_exp(0:jmaxmt), nm3_exp(0:jmaxmt)
      real*8 ptot_exp(0:jmaxmt), pfast_exp(0:jmaxmt)
      real*8 torque_exp(0:jmaxmt)
      real*8 vpar_exp(0:jmaxmt), vper_exp(0:jmaxmt)
      real*8 vphi_exp(0:jmaxmt), pzmn_sol(0:jmaxmt), ne_p(0:jmaxmt)
      real*8 vphi_ncl_exp(0:jmaxmt), vpol_ncl_exp(0:jmaxmt)
      real*8 omexb_ncl_exp(0:jmaxmt)
      real*8 ni_p(0:jmaxmt), zeff_exp(0:jmaxmt), zpte_exp(0:jmaxmt)
      real*8 zpti_exp(0:jmaxmt), zpne_exp(0:jmaxmt), zpni_exp(0:jmaxmt)
      real*8 zpnz_exp(0:jmaxmt), zpnitot_exp(0:jmaxmt)
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
      real*8 pow_ei_exp(0:jmaxmt), exch_exp(0:jmaxmt)
      real*8 flow_wall_exp(0:jmaxmt), flow_recom_exp(0:jmaxmt)
      real*8 flow_beam_exp(0:jmaxmt), flow_sdot_exp(0:jmaxmt)
      real*8 powineo_exp(0:jmaxmt)
      real*8 qbeame_exp(0:jmaxmt), qbeami_exp(0:jmaxmt)
      real*8 qrfe_exp(0:jmaxmt), qrfi_exp(0:jmaxmt)
      real*8 qohm_exp(0:jmaxmt), qrad_exp(0:jmaxmt)
      real*8 qione_exp(0:jmaxmt), qioni_exp(0:jmaxmt)
      real*8 sbeame_exp(0:jmaxmt), sbeam_exp(0:jmaxmt)
      real*8 flow_exch_exp(0:jmaxmt)
      real*8 powe_exp(0:jmaxmt), powi_exp(0:jmaxmt),flow_exp(0:jmaxmt)
      real*8 egamma_vphi(0:jmaxmt), egamma_vpol(0:jmaxmt)
      real*8 egamma_vstar(0:jmaxmt)
      real*8 chiegb_glf(0:jmaxmt), chiigb_glf(0:jmaxmt)
      real*8 diffgb_glf(0:jmaxmt), etagb_phi_glf(0:jmaxmt)
      real*8 etagb_par_glf(0:jmaxmt), etagb_per_glf(0:jmaxmt)
      real*8 rhosda_exp(0:jmaxmt)
      real*8 csda_exp(0:jmaxmt), bteff_exp(0:jmaxmt)
      real*8 vstar_exp(0:jmaxmt), vexb_exp(0:jmaxmt)
      real*8 vexbc_exp(0:jmaxmt), vexbcc_exp(0:jmaxmt)
      real*8 vstarp_exp(0:jmaxmt), egamma_exp(0:jmaxmt)
      real*8 egamma_ncl_exp(0:jmaxmt)
      real*8 megamma_exp(0:jmaxmt)
      real*8 gamma_p_i_exp(0:jmaxmt),gamma_p_e_exp(0:jmaxmt)
      real*8 mach_i_exp(0:jmaxmt),mach_e_exp(0:jmaxmt)
      real*8 gamma_mode_exp(0:jmaxmt), anrate_exp(0:jmaxmt)
      real*8 anfreq_exp(0:jmaxmt), cgyrobohm_exp(0:jmaxmt)
      real*8 betae_exp(0:jmaxmt), betai_exp(0:jmaxmt)
      real*8 alpha_exp(0:jmaxmt), xnu_exp(0:jmaxmt)
      real*8 vnewstare_exp(0:jmaxmt), vnewstari_exp(0:jmaxmt)
      real*8 chiegb_exp(0:jmaxmt), chiigb_exp(0:jmaxmt)
      real*8 diffgb_exp(0:jmaxmt), chiineogb_exp(0:jmaxmt)
      real*8 chiineo_exp(0:jmaxmt), chiicond_exp(0:jmaxmt)
      real*8 chieneo_exp(0:jmaxmt), chieneogb_exp(0:jmaxmt)
      real*8 chiecond_exp(0:jmaxmt)
      real*8 dineo_exp(0:jmaxmt), deneo_exp(0:jmaxmt)
      real*8 taue_exp(0:jmaxmt), taui_exp(0:jmaxmt), taut_exp(0:jmaxmt)
      real*8 taupe_exp(0:jmaxmt), taupi_exp(0:jmaxmt), rho(0:jmaxmt)
      real*8 rmin_exp(0:jmaxmt), rmaj_exp(0:jmaxmt), sfactor(0:jmaxmt)
      real*8 geofac(0:jmaxmt), drhodr(0:jmaxmt), drhodrrrho(0:jmaxmt)
      real*8 georotrate(0:jmaxmt), geoalpha(0:jmaxmt)
      real*8 h_exp(0:jmaxmt), f_exp(0:jmaxmt), g_exp(0:jmaxmt)
      real*8 xb2_exp(0:jmaxmt), xbm2_exp(0:jmaxmt), xngrth_exp(0:jmaxmt)
     & ,     xgrbm2_exp(0:jmaxmt), fm1_exp(0:jmaxmt), fm2_exp(0:jmaxmt)
     & ,     fm3_exp(0:jmaxmt), fhat_exp(0:jmaxmt)
      real*8 gradrhosq_exp(0:jmaxmt), gradrho_exp(0:jmaxmt)
      real*8 sfarea_exp(0:jmaxmt), cxarea_exp(0:jmaxmt)
      real*8 psir_exp(0:jmaxmt)
      real*8 vsurf_exp(0:jmaxmt)
      real*8 elong_exp(0:jmaxmt), delta_exp(0:jmaxmt)
      real*8 delt_exp(0:jmaxmt), q_exp(0:jmaxmt)
      real*8 curden_exp(0:jmaxmt), curboot_exp(0:jmaxmt)
      real*8 curbeam_exp(0:jmaxmt), currf_exp(0:jmaxmt)
      real*8 curohm_exp(0:jmaxmt)
      real*8 shat_exp(0:jmaxmt), er_exp(0:jmaxmt)
      real*8 angrot_exp(0:jmaxmt), angrotp_exp(0:jmaxmt)
      real*8 vecur_exp(0:jmaxmt)
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
c 100 elements
      real*8 te_cen_ei(1:icalleimx), te_mid_ei(1:icalleimx)
c
c 9 elements
      real*8 dflowdzp(3,3),rmatrix(3,3),one(3,3)
c 5 elements
      real*8 tauescale_exp(5)
c 3 elements
      real*8 zm(3), zpl(3), zph(3)
      real*8 zpsol(3),flowl(3),flowexp(3),fvector(3)
c
      real*8 kevdsecpmw
      real*8 pzmn(3),pzmn1,pzmn2
      real*8 detamatrix,relxei,relxbr
      real*8 relxfus,tau_he,f_defect,f_def_beryllium,relx
      real*8 eps_tetest,ascale,bscale,pscale,mscale,nfscale
      real*8 chisqte,chisqti,chisqne
      real*8 pbescale, pbiscale, snbscale
      real*8 prfscale, prfescale, prfiscale
      real*8 nscale
      real*8 sigmate,sigmati,sigmane,fte,fti,fne
      real*8 neline_exp,niline_exp
      real*8 zpmte1, zpmti1,zpmne1,zpmni1, volavene_exp
      real*8 volaveniti_exp,gtaut_exp,gtaui_exp
      real*8 gtauth_exp,gtautot_exp
      real*8 n_t_tau_exp,n0_t0_tau_exp,volaveniti2_exp
      real*8 wstre_exp,wstri_exp,wstre_ped_exp,wstri_ped_exp
      real*8 wstr_exp, wstr_inc_exp, wstr_core_exp, wstr_ped_exp
      real*8 p_glob_exp,w_glob_exp
      real*8 rhosda20_exp,te20_exp,ti20_exp,nebar_exp
      real*8 rho_g,btvac_toq,btplas_toq
      real*8 torfO2pi_toq,a_bnd_toq,a_unit_exp,bt_mag_center
      real*8 alpha_e_dk,alpha_star_dk
      real*8 xp_time, tol_f
      real*8 xfus_exp,xrad_exp,xoh_exp, xion_exp
      real*8 kys0,xalpha,alpha_p,alpha_p_cur,alpha_e,cbetae
      real*8 alpha_mach,alpha_cur
c
      integer nstk_s(150), mprint(10)
      integer iskip(stkmax), istk_smooth(stkmax)
      integer jstmix_s(stkmax)
      integer iei,iquit,iquitmax,ibifurcate,imodcall,jm,iloop_p_max
      integer iterate,iloopset,iiloopset,iiiloopset,iiiiloopset
      integer ivloopset,iscale,iexp_exch,iexp_q,iexp_neo
      integer iexp_imp,ineo,ineophi,idiffphi_dv
      integer jstmix,jstinv,ismooth,istep_smoo
      integer i_delay,i_gridsmooth,iv_dk,itest_dk, kmaxm
c
      character(50) cudir
      character(40) shot
      character(6) phase
      character(10) tok
      character(3) scanid
c
      common /data_exp/ sbcx_exp
     & , dineo_exp, deneo_exp
     & , te_exp, ti_exp, ne_exp, ni_exp, nz_exp, nfst_exp, nitot_exp
     & , nm1_exp, nm2_exp, nm3_exp
     & , ptot_exp, pfast_exp, vpar_exp, vper_exp 
     & , vphi_exp, vphi_ncl_exp, vpol_ncl_exp
     & , omexb_ncl_exp, pzmn_sol, ne_p
     & , ni_p, zeff_exp, zpte_exp, zpti_exp
     & , zpne_exp, zpni_exp, zpnitot_exp, zpti2_exp
     & , rlne_exp, rlni_exp, rlte_exp, rlti_exp
     & , zpnz_exp, vol_exp, powe_beam_exp, powe_rf_exp
     & , powe_lh_exp, powe_oh_exp
     & , powe_rad_exp, powe_ion_exp, powe_wdot_exp, powe_fus_exp
     & , powi_beam_exp, powi_rf_exp, powi_ion_exp, powi_cx_exp
     & , powi_wdot_exp, powi_fus_exp, pow_ei_exp, exch_exp
     & , flow_wall_exp, flow_recom_exp, flow_beam_exp
     & , flow_sdot_exp, powineo_exp
     & , flow_exch_exp, powe_exp, powi_exp, flow_exp
     & , egamma_vphi, egamma_vpol, egamma_vstar
     & , rhosda_exp
     & , csda_exp, bteff_exp, vstar_exp, vexb_exp, vexbc_exp
     & , vexbcc_exp, vstarp_exp, egamma_exp, egamma_ncl_exp
     & , megamma_exp, gamma_p_i_exp, gamma_p_e_exp
     & , mach_i_exp,mach_e_exp
     & , gamma_mode_exp, anrate_exp, anfreq_exp, cgyrobohm_exp
     & , betae_exp, betai_exp, alpha_exp, xnu_exp
     & , vnewstare_exp, vnewstari_exp, chiegb_exp, chiigb_exp
     & , diffgb_exp, chiineogb_exp, chiineo_exp
     & , chieneo_exp, chieneogb_exp
     & , chiicond_exp, chiecond_exp, taue_exp
     & , taui_exp, taut_exp, taupe_exp, taupi_exp, rho
     & , rmin_exp, rmaj_exp, sfactor, geofac, drhodr, drhodrrrho
     & , georotrate, geoalpha, h_exp, xb2_exp, xbm2_exp
     & , xngrth_exp, xgrbm2_exp,fm1_exp, fm2_exp, fm3_exp
     & , fhat_exp, gradrhosq_exp, gradrho_exp, vsurf_exp
     & , elong_exp, delta_exp, delt_exp, q_exp, curden_exp, shat_exp
     & , angrot_exp, angrotp_exp, vecur_exp,psi_exp, betat_exp
     & , xfus_exp,xrad_exp,xoh_exp, xion_exp
     & , b_unit_loc_s, b_norm_loc_s, rhosda_loc_s, csda_loc_s
     & , aspectratio_loc_s, kappa_loc_s, s_kappa_loc_s, delta_loc_s
     & , s_delta_loc_s, shift_loc_s, beta_loc_s, q_loc_s
     & , shat_loc_s, dlntidr_loc_s, dlntedr_loc_s, dlnnidr_loc_s
     & , dlnnedr_loc_s, tiote_loc_s, nione_loc_s, xnu_loc_s
     & , shat_mhd_loc_s, alpha_mhd_loc_s, alpha_mhd1, alpha_mhd2
     & , beta_loc_0_s, volume_loc_s, delta_dor_s, d_prime_dor_s
     & , k_prime_dor_s
     & , te_cen_ei, te_mid_ei
     & , dflowdzp, rmatrix, one, tauescale_exp
      common /data_exp/ zm, zpl, zph, zpsol,flowl,flowexp
     & , fvector, kevdsecpmw
     & , pzmn,pzmn1,pzmn2
     & , detamatrix,relxei,relxbr
     & , relxfus,tau_he,f_defect,f_def_beryllium,relx,eps_tetest
     & , ascale,bscale,pscale,mscale,nfscale
     & , pbescale,pbiscale,snbscale,prfscale,prfescale,prfiscale
     & , nscale,chisqte,chisqti
     & , chisqne,sigmate,sigmati,sigmane,fte,fti,fne
     & , neline_exp,niline_exp,zpmte1
     & , zpmti1,zpmne1,zpmni1,volavene_exp,volaveniti_exp
     & , gtaut_exp,gtaui_exp,gtauth_exp,gtautot_exp
     & , n_t_tau_exp,n0_t0_tau_exp,volaveniti2_exp,wstre_exp,wstri_exp
     & , wstre_ped_exp,wstri_ped_exp,p_glob_exp,w_glob_exp
     & , rhosda20_exp,te20_exp,ti20_exp,nebar_exp
     & , rho_g,btvac_toq,btplas_toq
     & , torfO2pi_toq,a_bnd_toq,a_unit_exp,bt_mag_center,alpha_e_dk
     & , alpha_star_dk, xp_time, tol_f
     & , wstr_exp, wstr_inc_exp, wstr_core_exp, wstr_ped_exp
     & , kys0,xalpha
     & , alpha_p,alpha_p_cur,alpha_e,cbetae,alpha_mach,alpha_cur
     & , nstk_s, mprint, iskip, istk_smooth, jstmix_s
     & , iei,iquit,iquitmax,ibifurcate,imodcall,jm,iloop_p_max
     & , iterate,iloopset,iiloopset,iiiloopset,iiiiloopset
     & , ivloopset,iscale,iexp_exch,iexp_q,iexp_neo
     & , iexp_imp,ineo,ineophi,idiffphi_dv
     & , jstmix,jstinv,ismooth,istep_smoo
     & , i_delay,i_gridsmooth,iv_dk,itest_dk,kmaxm
     & , cudir, shot, tok, phase, scanid

