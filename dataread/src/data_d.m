c@data_d.m 20-feb-02 J. Kinsey, General Atomics
c 25-feb-11 added gcap_d, hcap_d, psir_d, volsn_d, curboot_d,
c           curbeam_d, currf_d, rbp_d, bp0_d, cxareanpsi_d, er_d
c 20-feb-01 added ismoo_grho to smooth gradrho, gradrhosq
c 16-oct-01 added sbion_d for beam electron source
c 25-jan-01 added qlhe_d for lower hybrid
c 12-jan-01 added torque_d
c 14-nov-00 added en_nm1_d and en_nm2_d, en_nm3_d
c 16-may-00 added ipptot to include fast ions in ptot
c---------------------------------------------------------------------
c
      integer nplasb,nj,nion
      parameter (nplasb=4000, nj=301, nion=5)
c
      integer ishot_d, nj_d, nion_d, nprim_d, nimp_d,
     &   nneu_d, ibion_d, nplasbdry_d
c
      character(2) namep_d(nion), namei_d(nion), namen_d(nion)
c
      real*8 time_d, rgeom_d, rmag_d, rmajor_d, kappa_d,
     &   deltao_d, pindento_d, volo_d, areao_d, btor_d,
     &   tocur_d, totohm_d, totboot_d, totbeam_d, totrf_d,
     &   betap_d, beta_d, ali_d, te0_d, ti0_d,
     &   pohm_d, vsurf_d, amin_d, p_glob_d, gtauth_d, gtautot_d
c 1d arrays
      real*8 blank_d(nj), rho_d(nj), r_d(nj), rbp_d(nj), bp0_d(nj),
     &   fcap_d(nj), gcap_d(nj), hcap_d(nj), psir_d(nj),
     &   xb2_d(nj), xbm2_d(nj), xngrth_d(nj), xgrbm2_d(nj),
     &   fm1_d(nj), fm2_d(nj), fm3_d(nj), fhat_d(nj),
     &   te_d(nj), ti_d(nj), q_d(nj), ene_d(nj), enbeam_d(nj),
     &   sbion_d(nj), sbeam_d(nj), zeff_d(nj), angrot_d(nj),
     &   curden_d(nj), curohm_d(nj), curboot_d(nj), curbeam_d(nj),
     &   currf_d(nj),
     &   chieinv_d(nj), chiinv_d(nj), xkineo_d(nj),
     &   dpedtc_d(nj), dpidtc_d(nj), qconde_d(nj), 
     &   qcondi_d(nj), qconve_d(nj), qconvi_d(nj), 
     &   qbeame_d(nj), qdelt_d(nj), qbeami_d(nj),
     &   qrfe_d(nj), qrfi_d(nj), qlhe_d(nj),
     &   qione_d(nj), qioni_d(nj),
     &   qcx_d(nj), qfuse_d(nj), qfusi_d(nj), qrad_d(nj),
     &   qohm_d(nj), rmajavnpsi_d(nj), rminavnpsi_d(nj),
     &   psivolp_d(nj), elongx_d(nj), deltax_d(nj),
     &   sfareanpsi_d(nj), cxareanpsi_d(nj),
     &   grho1npsi_d(nj), grho2npsi_d(nj),
     &   er_d(nj), en_nm1_d(nj), en_nm2_d(nj), en_nm3_d(nj),
     &   torque_d(nj), ptot_d(nj), pfast_d(nj)
c.. 2d elements
      real*8 sion_d(nj,nion), srecom_d(nj,nion), scx_d(nj,nion),
     &   sbcx_d(nj,nion), s_d(nj,nion), dudtsv_d(nj,nion), 
     &   enn_d(nj,nion), ennw_d(nj,nion), ennv_d(nj,nion),
     &   volsn_d(nj,nion),en_d(nj,nion)
c.. 4000 elements
      real*8 bblank_d(nplasb)
c
      common /data_d0/
     &   ishot_d, nj_d, nion_d, nprim_d, nimp_d,
     &   nneu_d, ibion_d, nplasbdry_d,
     &   time_d, rgeom_d, rmag_d, rmajor_d, kappa_d,
     &   deltao_d, pindento_d, volo_d, areao_d, btor_d,
     &   tocur_d, totohm_d, totboot_d, totbeam_d, totrf_d,
     &   betap_d, beta_d, ali_d, te0_d, ti0_d,
     &   pohm_d, vsurf_d, amin_d, p_glob_d, gtauth_d, gtautot_d
      common /data_dc/ 
     &    namep_d, namei_d, namen_d
      common /data_d1/
     &   blank_d, rho_d, r_d, rbp_d, bp0_d,
     &   fcap_d, gcap_d, hcap_d, psir_d,
     &   xb2_d, xbm2_d, xngrth_d, xgrbm2_d,
     &   fm1_d, fm2_d, fm3_d, fhat_d,
     &   te_d, ti_d, q_d, ene_d, enbeam_d,
     &   sbion_d, sbeam_d, zeff_d, angrot_d,
     &   curden_d, curohm_d, curboot_d, curbeam_d,
     &   currf_d,
     &   chieinv_d, chiinv_d, xkineo_d,
     &   dpedtc_d, dpidtc_d, qconde_d, 
     &   qcondi_d, qconve_d, qconvi_d, 
     &   qbeame_d, qdelt_d, qbeami_d,
     &   qrfe_d, qrfi_d, qlhe_d,
     &   qione_d, qioni_d,
     &   qcx_d, qfuse_d, qfusi_d, qrad_d,
     &   qohm_d, rmajavnpsi_d, rminavnpsi_d,
     &   psivolp_d, elongx_d, deltax_d,
     &   sfareanpsi_d, cxareanpsi_d,
     &   grho1npsi_d, grho2npsi_d,
     &   er_d, en_nm1_d, en_nm2_d, en_nm3_d,
     &   torque_d, ptot_d, pfast_d
      common /data_d2/
     &   sion_d, srecom_d, scx_d,
     &   sbcx_d, s_d, dudtsv_d, 
     &   enn_d, ennw_d, ennv_d,
     &   volsn_d, en_d
      common /data_d3/ bblank_d

