module glf23_gf

  implicit none
  save

  logical :: use_transport_model_gf = .true.
  logical :: use_adiabatic_electrons_gf = .false.
  integer,parameter :: nsm=3

  integer :: version_gf=2  !re-tuned v1.61
  integer,parameter :: nmode=20
  integer,dimension(30) ::  iflagin_gf=0
  integer,dimension(0:nmode) :: ngrow_k_gf=0
  integer :: ns_gf = 2 
  integer :: nroot_gf=8
  integer :: nbasis_gf=4
  integer :: lprint_gf=0
  integer :: ikymax_gf=10
  integer :: eigen_gf=0
  integer :: i_err=0
  
  real :: yparam_k_gf(nmode,nmode)=0.0
  real :: gamma_k_gf(1:4,nmode)
  real :: freq_k_gf(1:4,nmode)
  real :: phi_norm_k_gf(1:4,nmode)
  real,dimension(30) :: xparam_gf=0.0
  real :: xkyf_k_gf(nmode)
  real :: diff_k_gf(nmode)
  real :: diff_im_k_gf(nmode)
  real :: chii_k_gf(nmode)
  real :: chi_im_k_gf(nmode)
  real :: chie_k_gf(nmode)
  real :: exch_k_gf(nmode)
  real :: eta_par_k_gf(nmode)
  real :: eta_per_k_gf(nmode)
  real :: eta_phi_k_gf(nmode)
  real :: chie_e_k_gf(nmode)
  real :: n_tilda2_k_gf(nmode)
  real :: yparam_gf(nmode)
  real :: gamma_gf(1:4)
  real :: freq_gf(1:4)
  real :: phi_norm_gf(1:4)
  real :: xky_gf(1:4)
  real :: particle_flux_gf(2,2)
  real :: energy_flux_gf(2,2)
  real :: xky0_gf = 0.3
  real :: rms_theta_gf=3.14149265/3.0
  real :: theta0_gf = 0.0
  real :: rlti_gf = 3.0
  real :: rlte_gf = 3.0
  real :: rlne_gf = 1.0
  real :: rlni_gf = 1.0
  real :: rlnimp_gf = 1.0
  real :: dil_gf = 0.0
  real :: apwt_gf = 1.0
  real :: aiwt_gf = 0.0
  real :: taui_gf = 1.0
  real :: xnu_gf = 0.0
  real :: betae_gf = 1.0E-12
  real :: rmin_gf = 0.5
  real :: rmaj_gf = 3.0
  real :: q_gf = 2.0
  real :: shat_gf = 1.0
  real :: alpha_gf = 0.0
  real :: elong_gf = 0.0
  real :: xwell_gf = 0.0
  real :: park_gf = 1.0
  real :: phii_gf
  real :: ghat_gf=1.0
  real :: gchat_gf=1.0
  real :: adamp_gf=0.7
  real :: alpha_star_gf = 0.0
  real :: gamma_star_gf = 0.0
  real :: alpha_e_gf = 1.0
  real :: alpha_e_mult_gf = 1.0
  real :: gamma_e_gf = 0.0
  real :: alpha_mode_gf = 0.0
  real :: gamma_mode_gf = 0.0
  real :: alpha_p_gf = 1.0
  real :: alpha_p_mult_gf=1.0
  real :: gamma_p_gf = 0.0
  real :: gamma_r_gf=0.0
  real :: xkdamp_gf=0.0
  real :: diff_gf
  real :: diff_im_gf
  real :: chii_gf
  real :: chi_im_gf
  real :: chie_gf
  real :: exch_gf
  real :: eta_par_gf
  real :: eta_per_gf
  real :: eta_phi_gf
  real :: chie_e_gf
  real :: cnorm_p_gf=1.0
  real :: cnorm_gf=1.0
  real :: xkymin_gf
  real :: xkymax_gf
  real :: amassgas_gf = 1.0
  real :: amassimp_gf = 6.0
  real :: zimp_gf = 6.0

  complex :: zevec_k_gf(nmode,12,12)
  complex :: zomega_k_gf(nmode,12)
  
end module glf23_gf


     
