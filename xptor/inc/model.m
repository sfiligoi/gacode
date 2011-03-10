c@model.m 17-Nov-05 J. Kinsey, General Atomics
c 17-nov-05 (jek) removed _sa,_gs variables
c 21-oct-03 added ntilda_m
c 09-oct-02 removed fgbe_m_s,fgbi_m_s,fbe_m_s,fbi_m_s
c          fge_m_s, fgi_m_s sized at (0:jmaxmm,stkmxm)
c 12-dec-01 added xion_exp for ionization losses
c 18-aug-01 added diff_k_j, chie_etg_j, chie_itg_j, etaphi_k_j
c 31-aug-00 added xchi,chiaddexp for artificial chi at small r
c 23-aug-00 added xkymax_m, xkymax2_m
c 11-aug-00 common block out of alignment, moved x_ff and f_ff
c 19-may-00 reorder common block due to increase in nj=51->301
c 1-05-00 added arrays for GLF23 gradients zpmte_glf, etc.
c 9-23-99 added cartdiff
c---------------------------------------------------------------------
      integer jmaxmm, stkmxm, nmxwn
      parameter (jmaxmm=300, stkmxm=150, nmxwn=9)
c..3010 elements
      real*8 gamma_k_j(10,0:jmaxmm), freq_k_j(10,0:jmaxmm)
      real*8 diff_k_j(10,0:jmaxmm), chie_k_j(10,0:jmaxmm)
      real*8 chii_k_j(10,0:jmaxmm), chie_itg_j(10,0:jmaxmm)
      real*8 chie_etg_j(10,0:jmaxmm), etaphi_k_j(10,0:jmaxmm)
c..1204 elements
      real*8 gamma_j(0:jmaxmm,1:4), freq_j(0:jmaxmm,1:4)
      real*8 phi_norm_j(0:jmaxmm,1:4)
c..450 elements
      real*8 x_ff(3,150),f_ff(3,150)
c..301 elements
      real*8 anrate_m(0:jmaxmm), anrate2_m(0:jmaxmm)
      real*8 anfreq_m(0:jmaxmm), anfreq2_m(0:jmaxmm)
      real*8 xkymax_m(0:jmaxmm), xkymax2_m(0:jmaxmm)
      real*8 dnrate_m(0:jmaxmm), dtnrate_m(0:jmaxmm)
      real*8 dnfreq_m(0:jmaxmm), zpti_dk(0:jmaxmm)
      real*8 zptim_dk(0:jmaxmm), zptim2_dk(0:jmaxmm)
      real*8 chi0_dk(0:jmaxmm), g0_dk(0:jmaxmm)
      real*8 ky_j(0:jmaxmm), zptineomax(0:jmaxmm)
      real*8 fgbe_m(0:jmaxmm), fgbi_m(0:jmaxmm)
      real*8 fbe_m(0:jmaxmm), fbi_m(0:jmaxmm)
      real*8 fge_m(0:jmaxmm), fgi_m(0:jmaxmm)
      real*8 fgbe_ave_m(0:jmaxmm), fgbi_ave_m(0:jmaxmm)
      real*8 fbe_ave_m(0:jmaxmm), fbi_ave_m(0:jmaxmm)
      real*8 fge_ave_m(0:jmaxmm), fgi_ave_m(0:jmaxmm)
      real*8 ntilda_m(0:jmaxmm)
      real*8 gne_wn(0:jmaxmm), gnh_wn(0:jmaxmm)
      real*8 gnz_wn(0:jmaxmm), gns_wn(0:jmaxmm)
      real*8 gte_wn(0:jmaxmm), gth_wn(0:jmaxmm)
      real*8 gtz_wn(0:jmaxmm), egamma_wn(0:jmaxmm)
      real*8 zpmte_glf(0:jmaxmm), zpmti_glf(0:jmaxmm)
      real*8 zpmne_glf(0:jmaxmm), zpmni_glf(0:jmaxmm)
c..81 elements
      real*8 difthi_wn(10,1:nmxwn)
c..32 elements
      real*8 cetain_wn(1:32)
c..10 elements
      real*8 amod0(10)
c..9 elements
      real*8 omega_wn(1:nmxwn)
      real*8 gamma_wn(1:nmxwn)
      real*8 velthi_wn(1:nmxwn)
      real*8 chieff_wn(1:nmxwn)
      real*8 chieff_wn17(1:nmxwn)
      real*8 perform_wn(1:nmxwn)
      real*8 flux_wn(1:nmxwn)
c..3 elements
      real*8 fig(3),fkb(3),frb(3)
c..1 element
      real*8 zpticrit,epsthres,thresitg0,dntem,dttem,ctteme
      real*8 cnteme,ctitge,cnitge,cttemi,cntemi,ctitgi,cnitgi,aconv
      real*8 etapar0,etaper0,etaphi0,xconv,eps_powtest,epsrate_m
      real*8 epsrate2_m,kys0,shift0,xalpha,alpha_p,cbetae,cxnu
      real*8 cmodel,cmodel_e,cmodel_i,cneo,cneophi,cneod
      real*8 xchi,chiaddexp,cartdiff
      real*8 alpha_neo,alpha_dia,a_rotstab
      real*8 xfus_m,xbr_m,xfus_exp,xrad_exp,xoh_exp, xion_exp
      real*8 dt_m,dtmin_m,ae_dk,rmajor_t
      real*8 epsilon_t,q_t,shat_t,betae_t,xnu_t,taui0_t,rlte_t,rlti_t
      real*8 rln_t,chiegb_t,chiigb_t,diffgb_t,anrate_t,dnrate_t
      real*8 dtnrate_t,anfreq_t,zeff_wn
      real*8 cetain_wn10,cetain_wn15,cetain_wn20,cetain_wn25
      real*8 cetain_wn27,cetain_wn30,cetain_wn32,sqrtelong_wn
      real*8 epsnhin_wn,epsnzin_wn,epsnsin_wn
      real*8 epstein_wn,epsthin_wn,epstzin_wn
      real*8 tauhin_wn,tauzin_wn,fnzin_wn
      real*8 czin_wn,azin_wn,fnsin_wn,betahin_wn,betazin_wn
      real*8 ftrapein_wn,ekyrhoin_wn,ekparlin_wn,zpmnh_wn, zpmni_wn
      real*8 zpmnz_wn,zpmns_wn,zpmte_wn,zpmth_wn,zpmtz_wn, zpmnq_wn
      real*8 nz_wn,tz_wn
      real*8 ky_wn,rhos_wn,omega_de_wn,betaein_wn,vefin_wn,qin_wn
      real*8 shearin_w,kappain_w,wexb_wn
      real*8 chie_wn_kb,chie_gd_kb,chie_sr_kb
      real*8 chii_wn_kb,chii_gd_kb,chii_sr_kb,diff_wn_kb,diff_gd_kb
      real*8 diff_sr_kb,zckb_kb,zcme_kb,zcmp_kb,zce_kb,vthe_kb,zgyrfe_kb
      real*8 zlare_kb,vthi_kb,zgyrfi_kb,zlari_kb,zlog_kb,zcf_kb,znuei_kb
      real*8 zep_kb,zlarpo_kb,zsglpr_kb,zlpr_kb,zsound_kb,zgdtime_kb
      real*8 zgdrlp_kb,zgdlna_kb,zrhos_kb,zgdln_kb,zgdalf_kb,zgddia_kb
      real*8 zgdp_kb,pow_tem_kb,zdtite_kb
      real*8 elong_kb,pow_elong_kb,cthery8_kb
      real*8 alphae_jt,alphai_jt,tol_f,stepmx_f
      real*8 eta_f,tol_n,stepmx_n,eps_n,consta_rl
      real*8 constc_rl,eta_rl,zeff_rl,curden_rl,critgrad_rl,gradte_rl
      real*8 gradti_rl,gradne_rl,gradq_rl,chian_rl,chii_rl,chie_rl
      real*8 heavyte_rl,heavyq_rl,rescale_rl,factr_rl
      real*8 mult_rlw
c
      integer jcount_f(0:jmaxmm)
      integer letain_wn(1:32)
      integer imodel,imodelfail,igbe,igbi,irot1,irot2,irotstab
      integer ialphastab,iexch_m,igeo_m,i_dengrad
      integer iexb,itest_exp
      integer iwn,iwnx,iwny,igridwn,itest_wn,nout_wn,letain_wn6
      integer letain_wn7,lprintin_wn,lprintin_gd,lprintin_kb,neq_wn
      integer ndim_wn,ngrad_wn,nmodes_wn,icond_wn,incl_wn_kb,incl_gd_kb
      integer incl_sr_kb,incl_elong_kb,jcount_f_ave,jcount_f_min
      integer jcount_f_max,count_f,maxcal_f,fitpowe_f,fitpowi_f
      integer fitflow_f,maxits_n,use_nag,use_xneo_m
      integer old_rlw,boucher_expr,add_cycl
      integer joop_flag
c
      common /modelcm/ gamma_k_j, freq_k_j
     & , diff_k_j, chie_k_j, chii_k_j
     & , chie_itg_j, chie_etg_j, etaphi_k_j
     & , gamma_j, freq_j, phi_norm_j, x_ff,f_ff
     & , anrate_m, anrate2_m, dnrate_m, dtnrate_m
     & , anfreq_m, anfreq2_m, xkymax_m, xkymax2_m
     & , dnfreq_m, zpti_dk, zptim_dk, zptim2_dk
     & , chi0_dk, g0_dk, ky_j, zptineomax
     & , fgbe_m, fgbi_m, fbe_m, fbi_m
     & , fge_m, fgi_m, fgbe_ave_m, fgbi_ave_m
     & , fbe_ave_m, fbi_ave_m, fge_ave_m, fgi_ave_m
     & , ntilda_m
     & , gne_wn, gnh_wn, gnz_wn, gns_wn
     & , gte_wn, gth_wn, gtz_wn, egamma_wn
     & , zpmte_glf, zpmti_glf, zpmne_glf, zpmni_glf
     & , difthi_wn, cetain_wn, amod0, omega_wn, gamma_wn
     & , velthi_wn, chieff_wn, chieff_wn17, perform_wn
     & , flux_wn, fig, fkb, frb
      common /modelcm/ zpticrit, epsthres, thresitg0
     & , dntem, dttem
     & , ctteme, cnteme, ctitge, cnitge, cttemi, cntemi
     & , ctitgi, cnitgi, aconv, etapar0, etaper0, etaphi0
     & , xconv, eps_powtest, epsrate_m, epsrate2_m, kys0
     & , shift0, xalpha, alpha_p, cbetae, cxnu, cmodel
     & , cmodel_e, cmodel_i, cneo, cneophi, cneod, xchi
     & , chiaddexp, cartdiff, alpha_neo
     & , alpha_dia, a_rotstab, xfus_m, xbr_m, xfus_exp
     & , xrad_exp, xoh_exp, xion_exp
     & , dt_m, dtmin_m, ae_dk, rmajor_t, epsilon_t
     & , q_t, shat_t, betae_t, xnu_t, taui0_t, rlte_t
     & , rlti_t, rln_t, chiegb_t, chiigb_t, diffgb_t
     & , anrate_t, dnrate_t, dtnrate_t, anfreq_t, zeff_wn
     & , cetain_wn10, cetain_wn15, cetain_wn20, cetain_wn25
     & , cetain_wn27, cetain_wn30, cetain_wn32, sqrtelong_wn
     & , epsnhin_wn, epsnzin_wn, epsnsin_wn, epstein_wn
     & , epsthin_wn, epstzin_wn, tauhin_wn, tauzin_wn
     & , fnzin_wn, czin_wn, azin_wn, fnsin_wn, betahin_wn
     & , betazin_wn, ftrapein_wn, ekyrhoin_wn, ekparlin_wn
     & , zpmnh_wn, zpmni_wn, zpmnz_wn, zpmns_wn
     & , zpmte_wn, zpmth_wn, zpmtz_wn, zpmnq_wn
     & , nz_wn, tz_wn, ky_wn
     & , rhos_wn, omega_de_wn, betaein_wn, vefin_wn, qin_wn
     & , shearin_w, kappain_w, wexb_wn, chie_wn_kb
     & , chie_gd_kb, chie_sr_kb, chii_wn_kb, chii_gd_kb
     & , chii_sr_kb, diff_wn_kb, diff_gd_kb, diff_sr_kb
     & , zckb_kb, zcme_kb, zcmp_kb, zce_kb, vthe_kb
     & , zgyrfe_kb, zlare_kb, vthi_kb, zgyrfi_kb, zlari_kb
     & , zlog_kb, zcf_kb, znuei_kb, zep_kb, zlarpo_kb
     & , zsglpr_kb, zlpr_kb, zsound_kb, zgdtime_kb
     & , zgdrlp_kb, zgdlna_kb, zrhos_kb, zgdln_kb, zgdalf_kb
     & , zgddia_kb, zgdp_kb, pow_tem_kb, zdtite_kb, elong_kb
     & , pow_elong_kb, cthery8_kb, alphae_jt, alphai_jt
     & , tol_f, stepmx_f, eta_f, tol_n, stepmx_n, eps_n
     & , consta_rl, constc_rl, eta_rl
     & , zeff_rl, curden_rl, critgrad_rl, gradte_rl
     & , gradti_rl, gradne_rl, gradq_rl, chian_rl, chii_rl
     & , chie_rl, heavyte_rl, heavyq_rl, rescale_rl
     & , factr_rl, mult_rlw
      common /modelcm/ jcount_f, letain_wn
     & , imodel, imodelfail, igbe, igbi, irot1, irot2
     & , irotstab, ialphastab, iexch_m, igeo_m, i_dengrad
     & , iexb, itest_exp, iwn, iwnx, iwny, igridwn
     & , itest_wn, nout_wn, letain_wn6, letain_wn7
     & , lprintin_wn, lprintin_gd, lprintin_kb, neq_wn
     & , ndim_wn, ngrad_wn, nmodes_wn, icond_wn, incl_wn_kb
     & , incl_gd_kb, incl_sr_kb, incl_elong_kb, jcount_f_ave
     & , jcount_f_min, jcount_f_max, count_f, maxcal_f
     & , fitpowe_f, fitpowi_f, fitflow_f, maxits_n, use_nag
     & , use_xneo_m, old_rlw, boucher_expr, add_cycl
     & , joop_flag	 

