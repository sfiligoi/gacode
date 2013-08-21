    MODULE tglfin
!-------------------------------------------------------------------------------
! copied from tglf_tgv193 
! It is extensively modified for use in gcnmp
!------------------------------------------------HSJ-2/05/09 update 6/30/2012---

      USE nrtype,                                                 ONLY : DP,I4B
!     USE tglf_dimensions,                                        ONLY : nsm
!      USE tglf_tgv193                          !tglfv1.93_modules.f90:
      USE tglf_tgv194                           !tglfv1.94_modules.f90


! ADDED by HSJ for compatibility:
!------------------------------------
 
     INTEGER,PARAMETER :: itprt_max = 5 
     INTEGER(I4B) tglf_max_pack         !HSJ
     INTEGER(I4B) ialphastab            !HSJ
     INTEGER(I4B) jt                    !HSJ use as an arbirary grid index for testing
    ! aliases for v1.93:
     LOGICAL      adi_elec_tg
     INTEGER(I4B) nb_max_tg,nb_min_tg
     INTEGER(I4B) nx_tg
     REAL(DP)     debye_fac_tg
     REAL(DP)     xky0_tg
     REAL(DP)     theta_trap_tg
     REAL(DP)     Linsker_tg
     REAL(DP)     gradB_tg
     REAL(DP)     etg_fac_tg
     REAL(DP)     tglf_version_no
     REAL(DP)     sign_It_tg_state     ! state file vlaues of signs +/- 1.
     REAL(DP)     sign_Bt_tg_state
    ! end aliases

     INTEGER(I4B)  nspecies_tg
     LOGICAL  tglf_unstable,monitor_stability
     LOGICAL       recall_tglf
     LOGICAL       test_tglf,tglf_startup_call,use_tglf_pert
     LOGICAL       tglf_grid_point_parallel      ! parallelize over grid
     LOGICAL       tglf_wave_no_parallel         ! parallelize over wave vector
     LOGICAL    tglf_implicit_grad
     REAL(DP)   tglf_flux_mult
     REAL(DP)   vpar_mach_tg
     REAL(DP)   vexb_mach_tg
     REAL(DP)   gamma_e_tg
     REAL(DP)   gamma_p_tg
     REAL(DP)   cnorm_tg
     REAL(DP)   xnuei_tg
     REAL(DP)   xnuei_fac_tg
     REAL(DP)   shift_tg
     REAL(DP)   x_psi_tg
     REAL(DP)   xalpha                  !HSJ
     REAL(DP)   xnu_tg
     REAL(DP)   rlne_tg
     REAL(DP)   rlni_tg
     REAL(DP)   rlnimp_tg
     REAL(DP)   rlte_tg
     REAL(DP)   rlti_tg
     REAL(DP)   aiwt_tg
     REAL(DP)   aimpwt_tg
     REAL(DP)   mass_ion_tg
     REAL(DP)   mass_imp_tg
     REAL(DP)   zion_tg
     REAL(DP)   zimp_tg
     REAL(DP)   xnu_fac_tg
     REAL(DP)   shat_mhd_tg
     REAL(DP)   alpha_mhd_tg
     REAL(DP)   amassgas_tg
     REAL(DP)   denfactor
     REAL(DP)   gb_unit
     REAL(DP)   alpha_dia   ! switch turns on off diamagnetic contribution to vexb
     REAL(DP)   zeffm
     REAL(DP)      tglf_rmin_sa_in,tglf_q_sa_in,tglf_shat_sa_in,                  &
                   tglf_alpha_sa_in,tglf_xwell_sa_in,tglf_theta0_sa_in


     ! for the multiple communicator group case we define some timing arrays
     REAL(DP),PUBLIC ::                                                            &
                                                 t_slave_end_tg,                   &
                                                 t_slave_start_tg,                 &
                                                 t_slave_compt_tg,                 &
                                                 t_start_tg,                       &
                                                 t_end_tg,                         &
                                                 t_master_com_tg

     INTEGER(I4B)  tglf_b_model_sa_in,tglf_ft_model_sa_in 

     CHARACTER*16, PARAMETER   ::  test_tglf_output ='test_tglf_output'
!--------------------------------------------------------------------------
! --  allocatable arrays overlap with  glf arrays but only one set is allocated
! --  in any given run. HSJ
! --  Note *_fg means full grid arrays, *_zc means zone centered arrays
!--------------------------------------------------------------------------
     REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:)             ::   &
          drho_drmin_fg,drho_drmin_zc,rmin_fg,rho_fg,rho_zc,   &
          drhodrrrho ,rmaj_fg,rmaj_zc,q_fg,q_zc,               & 
          zeff_zc,ne_fg ,ni_fg,nz_fg ,ns_fg,te_fg,ti_fg,       &
          vexb_zc,grad_ptot_zc,delta_fg,elong_fg,              &
          shat_zc,theta_fg,grad_vexb_zc,grad_nz_zc,            &
          vneo_zc,zptheta_zc,grad_vparm_zc,                    &
          fi_zc,fz_zc,betae_zc,betai_zc,beta_zc,               &
          drho_tg,alpha_zc,alpha_zc_save,                      &
          grad_ne_zc,grad_ni_zc,grad_te_zc,                    &
          egamma_zc,gamma_p_exp,dqdrho_zc,                     &
          egamma_zc_initial,csda_zc,grad_q_zc,                 &
          rhosda_fg, rhosda_zc,xnu_zc,cgyrobohm_fg,            &
          anrate_zc,anfreq_zc,csda_fg,b_unit_fg,               &
          nz_zc,ne_zc,ni_zc,te_zc,ti_zc,grad_ti_fg,            &
          grad_ion_press_fg,grad_ti_zc,tglf_send_packed,       &
          tglf_receive_packed
     REAL(DP),ALLOCATABLE,PUBLIC,SAVE,DIMENSION(:)        ::   & 
          tglf_typ_val 
 
                                        
     REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:,:,:)       ::   &
          tglf_p_flux,tglf_e_flux,tglf_m_flux

     REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:)           ::   &
          tglf_p_output,tglf_e_output,tglf_m_output,           &
          tglf_val_pert,tglf_val_base

     REAL(DP),ALLOCATABLE,PUBLIC,DIMENSION(:,:,:)         ::   &
          tglf_vp,tglf_vm,tglf_d
     REAL(DP) tglf_v_mult,tglf_d_mult

     INTEGER(I4B),ALLOCATABLE,PUBLIC,DIMENSION(:)         ::   &
          itport_pt
 
    END MODULE tglfin
