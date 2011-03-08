! tglfin.m   
! Created G. M. Staebler June 1, 2007
! fortran 77 style common block for TGLF input names read from the tglfin file
!
      INTEGER nrmax,nkmax,nxmax,nsmax
      PARAMETER(nrmax=300,nkmax=50,nxmax=47,nsmax=3)
!
! external switches
!
      INTEGER igeo_tg
      REAL*8 v2_bar
!
! single ky 
!
      REAL*8  xky0_tg
!
! eikonal 
!
      LOGICAL  new_eikonal_tg
!
! rare switches
!
      REAL*8  theta_trap_tg
      REAL*8  park_tg
      REAL*8  ghat_tg
      REAL*8  gchat_tg
      REAL*8  wd_zero_tg
      REAL*8  Linsker_tg
      REAL*8  gradB_tg
      REAL*8  filter_tg
      REAL*8  damp_psi_tg
      REAL*8  damp_sig_tg
!
! switches
!
      LOGICAL  iflux_tg
      LOGICAL  use_bper_tg
      LOGICAL  use_bpar_tg
      LOGICAL  use_mhd_rule_tg
      LOGICAL  use_bisection_tg
      INTEGER  ibranch_tg
      INTEGER  nmodes_tg
      INTEGER  nb_max_tg
      INTEGER  nb_min_tg
      INTEGER  nx_tg
      INTEGER  nky_tg
!
!  model parameters
!
      LOGICAL  adi_elec_tg
      REAL*8  alpha_p_tg
      REAL*8  alpha_e_tg
      REAL*8  alpha_kx0_tg
      REAL*8  alpha_kx1_tg
      REAL*8  alpha_kx2_tg
      REAL*8  alpha_quench_tg
      REAL*8  xnuei_fac_tg
      REAL*8  debye_fac_tg
      REAL*8  etg_fac_tg
      INTEGER sat_rule_tg
      INTEGER kygrid_model_tg
      INTEGER xnu_model_tg
      INTEGER vpar_model_tg      
!
! width 
!
      REAL*8  width_max_tg
      REAL*8  width_min_tg
      INTEGER  nwidth_tg
      LOGICAL  find_width_tg
!
! species
!
      INTEGER nspecies_tg
      REAL*8  zs_tg(nsmax)
      REAL*8  mass_tg(nsmax)
!
! gradients
!
      REAL*8  rlns_tg(nsmax)
      REAL*8  rlts_tg(nsmax)
      REAL*8  vpar_shear_tg(nsmax)
      REAL*8  vexb_shear_tg
!
! averages
!
      REAL*8  taus_tg(nsmax)
      REAL*8  as_tg(nsmax)
      REAL*8  vpar_tg(nsmax)
      REAL*8  betae_tg
      REAL*8  xnuei_tg
      REAL*8  zeff_tg
      REAL*8  debye_tg
!
!  common geometry
!
      REAL*8  rmin_tg
      REAL*8  rmaj_tg
      REAL*8  q_tg
!
! Shifted cicle inputs
!
      REAL*8  shat_tg
      REAL*8  alpha_tg
      REAL*8  xwell_tg
      REAL*8  theta0_tg
      INTEGER  b_model_tg
      INTEGER  ft_model_tg
!
! Miller inputs
!
      REAL*8  zmaj_tg
      REAL*8  drmindx_tg
      REAL*8  drmajdx_tg
      REAL*8  dzmajdx_tg
      REAL*8  kappa_tg
      REAL*8  s_kappa_tg
      REAL*8  delta_tg
      REAL*8  s_delta_tg
      REAL*8  zeta_tg
      REAL*8  s_zeta_tg
      REAL*8  q_prime_tg 
      REAL*8  p_prime_tg 
!
      COMMON /tglf_in/ adi_elec_tg,find_width_tg,new_eikonal_tg, 
     >   nb_max_tg,nb_min_tg,nx_tg,ibranch_tg,nspecies_tg, 
     >   nmodes_tg,iflux_tg,xky0_tg,width_max_tg,width_min_tg, 
     >   nwidth_tg,park_tg,ghat_tg,gchat_tg,alpha_e_tg,vexb_shear_tg, 
     >   alpha_p_tg,alpha_quench_tg,igeo_tg,theta_trap_tg, 
     >   theta0_tg,vpar_model_tg,v2_bar, 
     >   rmin_tg,rmaj_tg,zmaj_tg,use_bisection_tg, 
     >   q_tg,xnuei_tg,wd_zero_tg,betae_tg,shat_tg,alpha_tg, 
     >   xwell_tg,kappa_tg,s_kappa_tg,delta_tg,
     >   s_delta_tg,zeta_tg,s_zeta_tg, 
     >   drmindx_tg,drmajdx_tg,dzmajdx_tg,zeff_tg,debye_tg,use_bper_tg, 
     >   use_bpar_tg,use_mhd_rule_tg,q_prime_tg,damp_psi_tg,damp_sig_tg, 
     >   p_prime_tg,filter_tg,Linsker_tg,gradB_tg,  
     >   b_model_tg,ft_model_tg,xnuei_fac_tg,debye_fac_tg, 
     >   nky_tg,etg_fac_tg,kygrid_model_tg,xnu_model_tg, 
     >   sat_rule_tg,alpha_kx0_tg,alpha_kx1_tg,alpha_kx2_tg,   
     >   vpar_shear_tg,rlns_tg,rlts_tg,mass_tg,zs_tg,
     >   vpar_tg,taus_tg,as_tg,
