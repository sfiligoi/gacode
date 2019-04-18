module expro_locsim_interface

  integer :: n_species_exp
  
  double precision, parameter :: mass_deuterium  = 3.3452
  double precision, parameter :: temp_norm_fac   = 1602.2
  double precision, parameter :: charge_norm_fac = 1.6022

  double precision, dimension(:), allocatable :: rmin_exp

  double precision, dimension(:,:), allocatable :: temp_exp
  double precision, dimension(:,:), allocatable :: dens_exp
  double precision, dimension(:,:), allocatable :: dlntdr_exp
  double precision, dimension(:,:), allocatable :: dlnndr_exp
  double precision, dimension(:,:), allocatable :: sdlntdr_exp
  double precision, dimension(:,:), allocatable :: sdlnndr_exp

  double precision, dimension(:), allocatable :: gamma_e_exp
  double precision, dimension(:), allocatable :: gamma_p_exp
  double precision, dimension(:), allocatable :: mach_exp
  
  double precision, dimension(:,:,:), allocatable :: geo_yin_exp

  ! Local values

  double precision :: shift_loc
  double precision :: q_loc
  double precision :: s_loc
  double precision :: kappa_loc
  double precision :: delta_loc
  double precision :: zeta_loc
  double precision :: s_kappa_loc
  double precision :: s_delta_loc
  double precision :: s_zeta_loc
  double precision :: zmag_loc
  double precision :: dzmag_loc
  double precision :: gamma_e_loc
  double precision :: gamma_p_loc
  double precision :: mach_loc
  double precision :: rmin_loc
  double precision :: rmaj_loc
  double precision :: rhos_loc
  double precision :: z_eff_loc
  double precision :: b_unit_loc
  double precision :: rho_norm_loc
  double precision :: psi_norm_loc
  double precision :: psi_a_loc
  double precision :: cs_loc
  double precision :: betae_loc
  double precision :: beta_star_loc

  double precision, dimension(9) :: mass_loc
  double precision, dimension(9) :: z_loc
  double precision, dimension(9) :: dens_loc
  double precision, dimension(9) :: temp_loc
  double precision, dimension(9) :: dlnndr_loc
  double precision, dimension(9) :: dlntdr_loc
  double precision, dimension(9) :: sdlnndr_loc
  double precision, dimension(9) :: sdlntdr_loc

  integer :: geo_ny_loc
  double precision, dimension(:,:), allocatable :: geo_yin_loc

end module expro_locsim_interface
