subroutine expromake_init

  use EXPRO_interface

  ! Fundamental input.profiles scalars
  EXPRO_b_ref = 1.0
  EXPRO_arho  = 1.0

  ! Fundamental input.profiles arrays

  EXPRO_rho = 1.0
  EXPRO_rmin = 1.0
  EXPRO_rmaj = 3.0
  EXPRO_q = 2.0
  EXPRO_kappa = 1.0

  EXPRO_delta = 0.0
  EXPRO_te = 1.0
  EXPRO_ne = 1.0
  EXPRO_z_eff = 1.0
  EXPRO_w0 = 0.0

  EXPRO_flow_mom = 0.0
  EXPRO_pow_e    = 0.0
  EXPRO_pow_i    = 0.0
  EXPRO_pow_ei   = 0.0
  EXPRO_zeta     = 0.0

  EXPRO_flow_beam = 0.0
  EXPRO_flow_wall = 0.0
  EXPRO_zmag      = 0.0
  EXPRO_ptot      = 1.0
  EXPRO_poloidalfluxover2pi = 0.0

  EXPRO_ni = 1.0
  EXPRO_ti = 1.0
  EXPRO_vtor = 0.0
  EXPRO_vpol = 0.0

end subroutine expromake_init
