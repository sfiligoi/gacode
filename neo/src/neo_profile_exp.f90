module neo_profile_exp

  integer :: n_grid_exp
  integer :: n_species_exp

  real :: rhoN_torflux_a
  real, dimension(:), allocatable :: rhoN_torflux_exp
  real, dimension(:), allocatable :: rhoN_torflux
  real, dimension(:), allocatable :: psiN_polflux_exp
  real, dimension(:), allocatable :: rmin_exp
  real, dimension(:), allocatable :: rmaj_exp
  real, dimension(:), allocatable :: q_exp
  real, dimension(:), allocatable :: kappa_exp
  real, dimension(:), allocatable :: delta_exp
  real, dimension(:), allocatable :: zeta_exp
  real, dimension(:), allocatable :: zmag_exp

  real, dimension(:), allocatable :: te_ade_exp
  real, dimension(:), allocatable :: ne_ade_exp
  real, dimension(:), allocatable :: dlntdre_ade_exp
  real, dimension(:), allocatable :: dlnndre_ade_exp

  real, dimension(:,:), allocatable :: tem_exp
  real, dimension(:,:), allocatable :: den_exp

  real, dimension(:), allocatable :: r_p
  real, dimension(:), allocatable :: rmaj_p
  real, dimension(:), allocatable :: zmag_p
  real, dimension(:), allocatable :: b_unit_p
  real, dimension(:), allocatable :: shift_p
  real, dimension(:), allocatable :: s_kappa_p
  real, dimension(:), allocatable :: s_delta_p
  real, dimension(:), allocatable :: s_zeta_p
  real, dimension(:), allocatable :: s_zmag_p
  real, dimension(:), allocatable :: shat_p

  real, dimension(:,:), allocatable :: dlnndr_p
  real, dimension(:,:), allocatable :: dlntdr_p

  real, dimension(:), allocatable :: er_exp
  real, dimension(:), allocatable :: gradr0_p
  real, dimension(:), allocatable :: dphi0dr_p
  real, dimension(:), allocatable :: omega_rot_p
  real, dimension(:), allocatable :: omega_rot_deriv_p

  real, dimension(:,:,:), allocatable :: geo_yin_exp

end module neo_profile_exp


