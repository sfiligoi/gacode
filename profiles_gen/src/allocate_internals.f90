subroutine allocate_internals

  use prgen_globals

  implicit none

  allocate(rho(nx))
  allocate(dpsi(nx))
  allocate(rmin(nx))
  allocate(rmaj(nx))
  allocate(q(nx))
  allocate(p_tot(nx))
  allocate(omega0(nx))
  allocate(kappa(nx))
  allocate(delta(nx))
  allocate(zmag(nx))
  allocate(zeta(nx))
  allocate(vpolc_exp(nx))
  allocate(vtorc_exp(nx))

  allocate(johm(nx))
  allocate(jbs(nx))
  allocate(jnb(nx))
  allocate(jrf(nx))
 
end subroutine allocate_internals
