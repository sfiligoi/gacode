!---------------------------------------------------------------------
! Mapping from TGYRO internal variables to TGLF interface
!---------------------------------------------------------------------

subroutine tgyro_tglf_map

  use tgyro_globals
  use tglf_interface

  implicit none

  ! Local variables
  integer :: i_ion
  real :: q_abs

  ! Currently TGLF only works with positive q
  q_abs = abs(q(i_r))

  ! Initialize TGLF
  call tglf_init()

  !----------------------------------------------------------------
  ! Number of species (max=6)
  tglf_ns_in = loc_n_ion+1
  if (tglf_ns_in > nsm) then
     call tgyro_catch_error('ERROR: too many ions in TGLF')
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Charges: e,i,z
  tglf_zs_in(1) = -1.0
  do i_ion=1,loc_n_ion
     tglf_zs_in(i_ion+1) = zi_vec(i_ion)
  enddo
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Mass ratios: me/mi,mi/mi,mz/mi
  tglf_mass_in(1) = (me*loc_me_multiplier)/mi(1)
  do i_ion=1,loc_n_ion 
     tglf_mass_in(i_ion+1) = mi(i_ion)/mi(1)
  enddo
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! ne/ne,ni/ne,nz/ne
  tglf_as_in(1) = 1.0
  do i_ion=1,loc_n_ion
     tglf_as_in(i_ion+1) = ni(i_ion,i_r)/ne(i_r) 
  enddo
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Geometry parameters:
  !
  ! s-alpha
  tglf_rmin_sa_in     = r(i_r)/r_min
  tglf_rmaj_sa_in     = r_maj(i_r)/r_min
  tglf_q_sa_in        = q_abs
  tglf_shat_sa_in     = 1.0
  tglf_alpha_sa_in    = 0.0
  tglf_xwell_sa_in    = 0.0
  tglf_theta0_sa_in   = 0.0
  tglf_b_model_sa_in  = 0
  tglf_ft_model_sa_in = 1
  ! Miller
  tglf_rmin_loc_in    = r(i_r)/r_min
  tglf_rmaj_loc_in    = r_maj(i_r)/r_min
  tglf_q_loc_in       = q_abs
  tglf_q_prime_loc_in = (q_abs/(r(i_r)/r_min))**2*s(i_r)
  tglf_p_prime_loc_in = (q_abs/(r(i_r)/r_min))*(beta_unit(i_r)/(8*pi))*(-r_min*dlnpdr(i_r))
  tglf_shift_loc_in   = shift(i_r)
  tglf_kappa_loc_in   = kappa(i_r)
  tglf_s_kappa_loc_in = dmax1(s_kappa(i_r),1.e-2)
  tglf_delta_loc_in   = delta(i_r)
  tglf_s_delta_loc_in = s_delta(i_r)/sqrt(1.0-tglf_delta_loc_in**2)
  !----------------------------------------------------------------

  !-----------------------------------
  ! Density gradients (e,i,z)
  tglf_rlns_in(1) = r_min*dlnnedr(i_r)
  do i_ion=1,loc_n_ion
     tglf_rlns_in(i_ion+1) = r_min*dlnnidr(i_ion,i_r)
  enddo
  !-----------------------------------

  !-----------------------------------
  ! Temperature gradients (e,i,z)
  tglf_rlts_in(1) = r_min*dlntedr(i_r)
  do i_ion=1,loc_n_ion
     tglf_rlts_in(i_ion+1) = r_min*dlntidr(i_ion,i_r)
  enddo
  !-----------------------------------

  !-----------------------------------
  ! Te/Te,Ti/Te,Tz/Te
  tglf_taus_in(1) = 1.0
  do i_ion=1,loc_n_ion
     tglf_taus_in(i_ion+1) = ti(i_ion,i_r)/te(i_r)
  enddo
  !-----------------------------------

  !----------------------------------------------------------------
  ! Electron beta used for electromagnetic calculations
  tglf_betae_in = betae_unit(i_r)*loc_betae_scale
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Collisions:
  !
  ! Electron-ion collision frequency
  tglf_xnuei_in = nue(i_r)*r_min/c_s(i_r)*loc_nu_scale
  !
  ! Zeff
  tglf_zeff_in = z_eff(i_r)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Gamma_ExB (ExB shearing rate, units of a/cs)
  if (loc_rotation_method == 3) then
     if (loc_ebshear_flag == 1) tglf_vexb_shear_in = gamma_eb(i_r)*r_min/c_s(i_r)
     tglf_vpar_shear_in = gamma_p(i_r)*r_min/c_s(i_r)
  endif
  !----------------------------------------------------------------

  !-----------------------------------
  ! Number of high-k modes
  !  nky=12 (default to include ETG)
  !  nky=0  (low-k only)
  tglf_nky_in = 12
  !-----------------------------------

  !----------------------------------------------------------------
  ! Linear mode selection
  !
  ! i_branch_tg = -1: use nmodes_tg modes (default for transport)
  ! i_branch_tg =  0: most unstable mode
  ! i_branch_tg =  1: most unstable electron mode
  ! i_branch_tg =  2: most unstable ion mode
  tglf_ibranch_in = -1
  !
  ! Number of modes in transport calculation;
  ! should be 2-4.
  tglf_nmodes_in = 2
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Width selection parameters
  !
  ! Use bisection method to determine Gaussian width;
  ! bisection is better, faster.
  tglf_use_bisection_in = .true.
  !
  ! Number of Hermite functions to determine Gaussian Width (1-2)
  tglf_nbasis_min_in = 1
  !
  ! Number of Hermite function in full expansion (4 or more, even)
  tglf_nbasis_max_in = 4
  !
  ! Number of Hermite quadrature nodes
  tglf_nxgrid_in = 16
  !
  ! Maximum number of widths sampled (could be higher than 21)
  tglf_nwidth_in = 21
  !
  ! Bisection search interval; must increase nwidth_tg if
  ! increasing this interval to maintain accuracy:
  !
  !  accuracy ~ (width_max_tg-width_min_tg)/nwidth_tg
  !
  tglf_width_min_in = 0.3
  tglf_width_in = 1.65
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! CONTROL PARAMETERS
  !
  ! Include B_perp
  tglf_use_bper_in = .false.
  !
  ! Include B_parallel
  if (loc_betae_scale > 0.0) then
     tglf_use_bpar_in = .true.
  else
     tglf_use_bpar_in = .false.
  endif
  !
  ! Use adiabatic electrons
  tglf_adiabatic_elec_in = .false.
  !
  ! Compute eikonal (use stored value if false)
  tglf_new_eikonal_in = .true.
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! OTHER PARAMETERS
  !
  tglf_find_width_in    = .true.
  tglf_iflux_in         = .true.
  tglf_theta_trapped_in = 0.7
  tglf_xnu_factor_in    = 1.0
  tglf_debye_factor_in  = 1.0
  tglf_filter_in        = 2.0
  tglf_debye_in         = 0.0
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! New TGLF settings
  !
  select case (tgyro_tglf_revision)

  case (1)

     ! APS07 

     tglf_alpha_quench_in = 1.0
     tglf_alpha_e_in      = 0.0
     tglf_alpha_p_in      = 0.0
     tglf_alpha_kx0_in    = 0.0
     tglf_sat_rule_in     = 0
     tglf_kygrid_model_in = 1
     tglf_xnu_model_in    = 1

  case (2)

     ! Summer 2009

     tglf_alpha_quench_in = 1.0
     tglf_alpha_e_in      = 0.0
     tglf_alpha_p_in      = 0.0
     tglf_alpha_kx0_in    = 0.0
     tglf_sat_rule_in     = 0
     tglf_kygrid_model_in = 1
     tglf_xnu_model_in    = 2

  case (3)

     ! TGLF-09 w/o ExB shear

     tglf_alpha_quench_in = 0.0
     tglf_alpha_e_in      = 0.0
     tglf_alpha_p_in      = 0.0
     tglf_alpha_kx0_in    = 0.0
     tglf_sat_rule_in     = 0
     tglf_kygrid_model_in = 1
     tglf_xnu_model_in    = 2

  end select

  !----------------------------------------------------------------

  ! Dump parameters
  if (tgyro_tglf_dump_flag == 0) then
     tglf_dump_flag_in   = .false.
  else
     tglf_dump_flag_in = .true.
  endif
  tglf_dump_suffix_in = achar(i_proc_global+iachar("0"))
  tglf_dump_path_in   = ''

end subroutine tgyro_tglf_map
