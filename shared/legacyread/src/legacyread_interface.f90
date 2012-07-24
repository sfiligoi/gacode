!@data_interface.f90 March 18, 2011 G. Staebler, General Atomics
!---------------------------------------------------------------------

module data_interface

  implicit none

  integer, parameter :: nplasb=4000, nj=301, nion=5

  integer :: ishot_d, nj_d, nion_d, nprim_d, nimp_d,  &
       nneu_d, ibion_d, nplasbdry_d

  character(2) :: namep_d(nion), namei_d(nion), namen_d(nion)

  real :: time_d, rgeom_d, rmag_d, rmajor_d, kappa_d,  &
       deltao_d, pindento_d, volo_d, areao_d, btor_d,    &
       tocur_d, totohm_d, totboot_d, totbeam_d, totrf_d, &
       betap_d, beta_d, ali_d, te0_d, ti0_d,             &
       pohm_d, vsurf_d, amin_d, p_glob_d, gtauth_d, gtautot_d

  ! 1d arrays
  real :: blank_d(nj), rho_d(nj), r_d(nj), rbp_d(nj), bp0_d(nj), &
       fcap_d(nj), gcap_d(nj), hcap_d(nj), psir_d(nj),       &
       xb2_d(nj), xbm2_d(nj), xngrth_d(nj), xgrbm2_d(nj),    &
       fm1_d(nj), fm2_d(nj), fm3_d(nj), fhat_d(nj),      &
       te_d(nj), ti_d(nj), q_d(nj), ene_d(nj), enbeam_d(nj),    &
       sbion_d(nj), sbeam_d(nj), zeff_d(nj), angrot_d(nj),    &
       curden_d(nj), curohm_d(nj), curboot_d(nj), curbeam_d(nj),  &
       currf_d(nj), chieinv_d(nj), chiinv_d(nj), xkineo_d(nj),   &
       dpedtc_d(nj), dpidtc_d(nj), qconde_d(nj),     &
       qcondi_d(nj), qconve_d(nj), qconvi_d(nj),     &
       qbeame_d(nj), qdelt_d(nj), qbeami_d(nj),      &
       qrfe_d(nj), qrfi_d(nj), qlhe_d(nj), qione_d(nj), qioni_d(nj), &
       qcx_d(nj), qfuse_d(nj), qfusi_d(nj), qrad_d(nj),  &
       qohm_d(nj), rmajavnpsi_d(nj), rminavnpsi_d(nj),   &
       psivolp_d(nj), elongx_d(nj), deltax_d(nj),   &
       sfareanpsi_d(nj), cxareanpsi_d(nj),   &
       grho1npsi_d(nj), grho2npsi_d(nj),   &
       er_d(nj), en_nm1_d(nj), en_nm2_d(nj), en_nm3_d(nj),  &
       torque_d(nj), ptot_d(nj), pfast_d(nj)  

  !.. 2d elements
  real :: sion_d(nj,nion), srecom_d(nj,nion), scx_d(nj,nion), &
       sbcx_d(nj,nion), s_d(nj,nion), dudtsv_d(nj,nion),   &
       enn_d(nj,nion), ennw_d(nj,nion), ennv_d(nj,nion),   &
       volsn_d(nj,nion),en_d(nj,nion)

  !.. 4000 elements
  real :: bblank_d(nplasb)

end module data_interface
