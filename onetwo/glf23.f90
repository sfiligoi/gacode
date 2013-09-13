      Module glf23
!
!
!    jmaxm is used both as a dimension in glf23 and as the size
!    of the working arrays. For 1d arrays this is not a
!    problem but for the single 2d array, egamma_d, that is passed to
!    glf this requires special treatment. Hence this include file was turned
!    into a module instead. HSJ 03/25/02
       USE param,               ONLY : kion
       logical dv_method
       integer,parameter :: itport_glf  = 5
       integer ibtflag_glf,include_glf,limp_glf,ivphi_glf,     &
               lprint_glf,       &
              glf23_ncpus,jroot_iglf,glf23_iglf,               &
              jeigen_iglf,irotstab,i_delay,write_glf_namelist, &
              itte_dv,itti_dv,itangrot_dv,itenp_dv,            &
              init_egamma_exp,itene_dv
       integer, save, dimension(:), allocatable:: itport_pt
       real *8  exbmult_glf,x_alpha_glf,dv_delt
       real *8 iglf_chie,iglf_chii,iglf_chiv,iglf_d
       real *8 chid_glf_d,chid_glf_e,ichid_smooth !jmp.den
       real *8 chiv_glf_d,chiv_glf_e !jmp.den
       integer ichiv_model,ichiv_chii !jmp.snu
       integer iglf_eq !jmp.den
       real *8 iglf_idt !jmp.den
       real *8, save,dimension(:,:), allocatable:: egamma_d      
       real *8, save,dimension(:), allocatable::               &
              te_glf,ti_glf,ne_glf, ni_glf,                    &
              q_glf,  shat_exp,  elong_glf,  zeff_exp,         &
              rho_glf,  grho1_glf,  grho2_glf,  gradte_glf,    &
              zpte_glf,  gradti_glf,  zpti_glf,                &
              gradne_glf,  zpne_glf,  gradni_glf,  tau_glf,    &
              zpni_glf,  betae_glf,  betai_glf,  beta_glf,     &
              alpha_exp,  rmin_glf,  rmaj_glf,  drhogf,        &
              ns_glf,  bteff_glf,                              &
              csda_glf,  rhosda_glf,  drhodr,   drhodrrrho,    &
              angrotp_exp,  egamma_exp,  gamma_p_exp,          &
              vphim_glf,  vparm_glf,  vperm_glf,  ve_glf,      &
              fcglf,  ak1glf,  ak2glf,  aneoglf,               &
              diff_m,  chie_m,  chii_m,chie_etg_m,             &
              etaphi_m,  etapar_m,  etaper_m,                  &
              exch_m , egamma_m ,                              &
               gamma_p_m ,anrate_m , anrate2_m ,               &
              anfreq_m , anfreq2_m , gamma_net_i,              &
              gamma_net_e,                                     &
              ve_eff_m,vi_eff_m,vrot_eff_m,                    &
              diff_mdv,  chie_mdv,  chii_mdv,                  &
              etaphi_mdv,  etapar_mdv,  etaper_mdv,            &
              exch_mdv , egamma_mdv,egamma_exp_initial,        &
               gamma_p_mdv ,anrate_mdv , anrate2_mdv ,         &
              anfreq_mdv , anfreq2_mdv,zpti_save,              &
              chii_save,vexb_den_grad,vexb_ti_grad,            &
              vexb_w_grad,gamma_den,gamma_ti,gamma_w,          &
              vexb_tot
!
      real *8 chie_dv(-1:1),chii_dv(-1:1),chitorot_dv(-1:1),   &
              zptij_p(-1:1)
 
 
        END MODULE  glf23
