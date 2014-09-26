      MODULE glf23_gf
!
! replaces common glf.m
!
!---------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
!
      LOGICAL,EXTERNAL :: glf23_isnan
      LOGICAL,EXTERNAL :: glf23_isinf
      LOGICAL :: use_transport_model_gf = .true.
      LOGICAL :: use_adiabatic_electrons_gf = .false.
      INTEGER,PARAMETER :: nsm=3
!
      INTEGER :: version_gf=2  !re-tuned v1.61
      INTEGER,PARAMETER :: nmode=20
      INTEGER,DIMENSION(30) ::  iflagin_gf=0
      INTEGER,DIMENSION(0:nmode) :: ngrow_k_gf=0
      INTEGER :: ns_gf = 2 
      INTEGER :: nroot_gf=8
      INTEGER :: nbasis_gf=4
!      INTEGER :: nx_gf
!      INTEGER :: jeigen=0
!      INTEGER :: igfac
      INTEGER :: lprint_gf=0
      INTEGER :: ikymax_gf=10
      INTEGER :: eigen_gf=0
      INTEGER :: i_err=0
!      INTEGER :: first_order_gf
!      INTEGER :: ipert_gf=0
      INTEGER :: iglf=98
!      INTEGER :: ipade_gf=0
!      INTEGER :: ibranch_gf=0
!      INTEGER :: igauss_gf
!      INTEGER :: iflux_gf
!      INTEGER :: igeo_gf
!      INTEGER :: iroot_ion(2)
!      INTEGER :: iroot_elec(2)
!
      REAL :: yparam_k_gf(nmode,nmode)=0.0
      REAL :: gamma_k_gf(1:4,nmode)
      REAL :: freq_k_gf(1:4,nmode)
      REAL :: phi_norm_k_gf(1:4,nmode)
      REAL,DIMENSION(30) :: xparam_gf=0.0
      REAL :: xkyf_k_gf(nmode)
      REAL :: diff_k_gf(nmode)
      REAL :: diff_im_k_gf(nmode)
      REAL :: chii_k_gf(nmode)
      REAL :: chie_k_gf(nmode)
      REAL :: exch_k_gf(nmode)
      REAL :: eta_par_k_gf(nmode)
      REAL :: eta_per_k_gf(nmode)
      REAL :: eta_phi_k_gf(nmode)
      REAL :: chie_e_k_gf(nmode)
      REAL :: n_tilda2_k_gf(nmode)
      REAL :: yparam_gf(nmode)
      REAL :: gamma_gf(1:4)
      REAL :: freq_gf(1:4)
      REAL :: phi_norm_gf(1:4)
      REAL :: xky_gf(1:4)
      REAL :: particle_flux_gf(2,2)
      REAL :: energy_flux_gf(2,2)
      REAL :: xky0_gf = 0.3
      REAL :: rms_theta_gf=3.14149265/3.0
      REAL :: theta0_gf = 0.0
      REAL :: rlti_gf = 3.0
      REAL :: rlte_gf = 3.0
      REAL :: rlne_gf = 1.0
      REAL :: rlni_gf = 1.0
      REAL :: rlnimp_gf = 1.0
      REAL :: dil_gf = 0.0
      REAL :: apwt_gf = 1.0
      REAL :: aiwt_gf = 0.0
      REAL :: taui_gf = 1.0
      REAL :: xnu_gf = 0.0
      REAL :: betae_gf = 1.0E-12
      REAL :: rmin_gf = 0.5
      REAL :: rmaj_gf = 3.0
      REAL :: q_gf = 2.0
      REAL :: shat_gf = 1.0
      REAL :: alpha_gf = 0.0
      REAL :: elong_gf = 0.0
      REAL :: xwell_gf = 0.0
      REAL :: park_gf = 1.0
      REAL :: phii_gf
      REAL :: ghat_gf=1.0
      REAL :: gchat_gf=1.0
      REAL :: adamp_gf=0.7
      REAL :: alpha_star_gf = 0.0
      REAL :: gamma_star_gf = 0.0
      REAL :: alpha_e_gf = 1.0
      REAL :: alpha_e_mult_gf = 1.0
      REAL :: gamma_e_gf = 0.0
      REAL :: alpha_mode_gf = 0.0
      REAL :: gamma_mode_gf = 0.0
      REAL :: alpha_p_gf = 1.0
      REAL :: alpha_p_mult_gf=1.0
      REAL :: gamma_p_gf = 0.0
      REAL :: gamma_r_gf=0.0
      REAL :: xkdamp_gf=0.0
!      REAL :: xkyf_gf = 0.3
      REAL :: diff_gf
      REAL :: diff_im_gf
      REAL :: chii_gf
      REAL :: chie_gf
!      REAL :: diff_ib_gf
!      REAL :: diff_eb_gf
!      REAL :: chii_ib_gf
!      REAL :: chii_eb_gf
!      REAL :: chie_ib_gf
!      REAL :: chie_eb_gf
      REAL :: exch_gf
      REAL :: eta_par_gf
      REAL :: eta_per_gf
      REAL :: eta_phi_gf
      REAL :: chie_e_gf
!      REAL :: n_tilda_gf
      REAL :: cnorm_p_gf=1.0
      REAL :: cnorm_gf=1.0
!      REAL :: cnorme_gf
!      REAL :: cnormi_gf
      REAL :: xkymin_gf
      REAL :: xkymax_gf
      REAL :: amassgas_gf = 1.0
      REAL :: amassimp_gf = 6.0
      REAL :: zimp_gf = 6.0
!      REAL :: cnormi_tg
!      REAL :: cnorme_tg
!      REAL :: cnormd_tg
!      REAL :: chii_ib_tg
!      REAL :: chie_ib_tg
!      REAL :: chii_eb_tg
!      REAL :: chie_eb_tg
!      REAL :: diff_ib_tg
!      REAL :: diff_eb_tg

      COMPLEX :: zevec_k_gf(nmode,12,12)
      COMPLEX :: zomega_k_gf(nmode,12)
!
     END MODULE glf23_gf


     
