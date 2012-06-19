c@data.m J. Kinsey, General Atomics
c 25-aug-09 added volsn_d
c 19-aug-08 added cxareanpsi_d
c 21-feb-08 added er_d
c 20-feb-07 added pfast_d, ptot_d
c 16-dec-05 isolated neoclassical variables
c 20-feb-01 added ismoo_grho to smooth gradrho, gradrhosq
c 16-oct-01 added sbion_d for beam electron source
c 25-jan-01 added qlhe_d for lower hybrid
c 12-jan-01 added torque_d
c 14-nov-00 added en_nm1_d and en_nm2_d, en_nm3_d
c 13-nov-00 added nr_xp, rho_ncl_r to save #pts and rho from EFIT in FORCEBAL
c           also added itest_ncl to print neocl. transport from NCLASS via readxp.f
c 06-nov-00 added apgasa, apgasz, abgasa, abgasz, apimpa, apimpz
c 06-nov-00 added itest_ntcc for NTCC benchmarking (use kappa from 1D ufile)
c 08-aug-00 added ECH conversion factor, echconv
c 17-aug-00 added smoothing variables lprint_smoo, ismoo_ne, etc.
c 03-aug-00 added _new variables used to interpolate _d variables
c           in readxp.f when grid > 50. Also, declared rho_d
c 19-may-00 reorder common block due to increase in nj=51->301
c 16-may-00 added ipptot to include fast ions in ptot
c 03-may-00 added pimpa and pimpz variables
c 29-sept-99 added NCLASS variables
c---------------------------------------------------------------------
c
      character*(2) namep_d(3), namei_d(3), namen_d(3)
      integer ishot_d, nj_d, nj, nion_d, nprim_d, nimp_d,
     &   nneu_d, ibion_d, nplasbdry_d, nr_xp, ismooth_all, shotno,
     &   ipptot, lprint_smoo
      integer nplasb, ng_nc,
     &   numzones_nc, nkimod_nc, istringer_nc, ng,
     &   imyneoclass, itest_ntcc, itest_ncl
c
      parameter (nj=301, ng=3)
c
      real*8 time_d, rgeom_d, rmag_d, rmajor_d, kappa_d,
     &   deltao_d, pindento_d, volo_d, areao_d, btor_d,
     &   tocur_d, totohm_d, totboot_d, totbeam_d, totrf_d,
     &   betap_d, beta_d, ali_d, te0_d, ti0_d,
     &   pohm_d, pfusion_m, pbr_tot, vsurf_d, amin_d
      real*8 pimpa, pimpz,
     &   apgasa, apgasz, abgasa, abgasz, apimpa, apimpz
      real*8 ismoo_ne, ismoo_ni, ismoo_nf, ismoo_zeff, ismoo_q,
     &   ismoo_vrot, ismoo_qnb, ismoo_torque, ismoo_grho,
     &   ismoo_delta
      real*8 echconv
c
c... neoclassical variables
c
      real*8 btf_nc, drshaf_nc, rminor_nc, rmajor_nc,
     &   aimp_nco, xzimp_nco
      real*8 rhoel_nco, rhi_nco, rhoi_nco,
     &   te_nco, ti_nco, zeff_nco, q_nco, xkapi_nco,
     &   xnstari_nco, xnstare_nco, zfluxlim_nco
      real*8 aplasm_nc(ng), rhob_nco(ng)
      real*8 rho_ncl_r(130)
      real*8 rhoel_nc(nj), rhi_nc(nj), rhoi_nc(nj),
     &   te_nc(nj), ti_nc(nj), zeff_nc(nj), q_nc(nj),
     &   aimp_nc(nj), xzimp_nc(nj)
c
      real*8 rho_d(nj), r_d(nj), rbp_d(nj), bp0_d(nj),
     &   hcap_d(nj), fcap_d(nj), gcap_d(nj), xr2_d(nj), psir_d(nj),
     &   xb2_d(nj), xbm2_d(nj), xngrth_d(nj), xgrbm2_d(nj),
     &   fm1_d(nj), fm2_d(nj), fm3_d(nj), fhat_d(nj),
     &   te_d(nj), ti_d(nj), q_d(nj), ene_d(nj), enbeam_d(nj),
     &   pfast_d(nj), ptot_d(nj), zeff_d(nj), 
     &   sbion_d(nj), sbeam_d(nj),
     &   angrot_d(nj), curden_d(nj), curohm_d(nj), 
     &   curboot_d(nj), curbeam_d(nj), currf_d(nj),
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
     &   er_d(nj)
c
      real*8 en_nm1_d(nj), en_nm2_d(nj), en_nm3_d(nj),
     &   torque_d(nj)
      real*8 sion_d(nj,ng), srecom_d(nj,ng), scx_d(nj,ng),
     &   sbcx_d(nj,ng), s_d(nj,ng), dudtsv_d(nj,ng), 
     &   enn_d(nj,ng), ennw_d(nj,ng), ennv_d(nj,ng), volsn_d(nj,ng)
c.. 903 elements
c     real*8 rhob_nc(nj,ng)
c.. 1204 elements
      real*8 en_d(nj,ng+1)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      common /data/ en_d, sion_d, srecom_d, scx_d,
     &   sbcx_d, s_d, dudtsv_d, enn_d, ennw_d, ennv_d, volsn_d,
     &   rho_d, r_d, rbp_d, bp0_d, hcap_d, fcap_d, gcap_d, 
     &   xr2_d, psir_d, xb2_d, xbm2_d,
     &   xngrth_d, xgrbm2_d, fm1_d, fm2_d, fm3_d, fhat_d,
     &   te_d, ti_d, q_d, ene_d, enbeam_d, pfast_d, ptot_d,
     &   zeff_d, sbion_d, sbeam_d, curden_d, curohm_d, 
     &   curboot_d, curbeam_d, currf_d, angrot_d,
     &   chieinv_d, chiinv_d, xkineo_d, dpedtc_d, dpidtc_d,
     &   qconde_d, qcondi_d, qconve_d, qconvi_d, qbeame_d,
     &   qdelt_d, qbeami_d, qrfe_d, qrfi_d, qlhe_d,
     &   qione_d, qioni_d, qcx_d, qfuse_d, qfusi_d, qrad_d,
     &   qohm_d, rmajavnpsi_d, rminavnpsi_d, psivolp_d,
     &   elongx_d, deltax_d, sfareanpsi_d, cxareanpsi_d, 
     &   grho1npsi_d, grho2npsi_d, er_d
      common /data/ en_nm1_d, en_nm2_d, en_nm3_d, torque_d, 
     &   time_d, rgeom_d, rmag_d, rmajor_d, amin_d,
     &   kappa_d, deltao_d, pindento_d, volo_d, areao_d,
     &   btor_d, tocur_d, totohm_d, totboot_d, totbeam_d, 
     &   totrf_d, betap_d, beta_d, ali_d, te0_d, ti0_d,
     &   pohm_d, pfusion_m, pbr_tot, vsurf_d,
     &   ismoo_ne, ismoo_ni, ismoo_nf, ismoo_zeff, ismoo_q,
     &   ismoo_vrot, ismoo_qnb, ismoo_torque, ismoo_grho,
     &   ismoo_delta, echconv, pimpa, pimpz,
     &   apgasa, apgasz, abgasa, abgasz, apimpa, apimpz,
     &   ishot_d, nj_d, nion_d, nprim_d, nimp_d, ismooth_all,
     &   shotno, nneu_d, ibion_d, nplasbdry_d, 
     &   imyneoclass, ipptot,
     &   lprint_smoo, itest_ntcc, itest_ncl, nr_xp,
     &   namep_d, namei_d, namen_d
c
      common /dataneo/ rhoel_nc, rhi_nc, rhoi_nc,
     &   te_nc, ti_nc, zeff_nc, q_nc, aimp_nc, xzimp_nc, rho_ncl_r,
     &   aplasm_nc, rhob_nco, 
     &   btf_nc, drshaf_nc, rminor_nc, rmajor_nc,
     &   aimp_nco, xzimp_nco,
     &   rhoel_nco, rhi_nco, rhoi_nco,
     &   te_nco, ti_nco, zeff_nco, q_nco, xkapi_nco,
     &   xnstari_nco, xnstare_nco, zfluxlim_nco,
     &   ng_nc, numzones_nc, nkimod_nc, istringer_nc
