subroutine expromake_init

  use expromake_globals
  use EXPRO_interface

  implicit none

  integer :: i
  real :: x


  ! Fundamental input.profiles scalars
  EXPRO_b_ref = exm_b_ref
  EXPRO_arho  = exm_arho

  ! Fundamental input.profiles arrays

  do i=1,nx
     x = real(i-1)/(nx-1) 
     EXPRO_rho(i)  = x
     EXPRO_rmin(i) = 0.6*x 
     EXPRO_te(i)   = exm_te_axis*exp(-x*exm_alte)
  enddo

  EXPRO_rmaj = 1.7
  EXPRO_q = 2.0
  EXPRO_kappa = exm_kappa
  EXPRO_delta = exm_delta

  EXPRO_ne = 5.0
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

  EXPRO_ni = 3.0
  EXPRO_ti = 7.0
  EXPRO_vtor = 0.0
  EXPRO_vpol = 0.0

end subroutine expromake_init
