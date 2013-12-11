subroutine allocate_internals

  use prgen_globals

  implicit none


  ! Profiles
 
  allocate(rmin(nx))
  allocate(rmaj(nx))
  allocate(q(nx))
  allocate(rho(nx))
  allocate(dpsi(nx))
  allocate(kappa(nx))
  allocate(delta(nx))
  allocate(zmag(nx))
  allocate(zeta(nx))
  allocate(omega0(nx))
  allocate(q_gato(nx))

   ! Powers

  allocate(powe_beam_exp(nx))
  allocate(powi_beam_exp(nx))
  allocate(powe_rf_exp(nx))
  allocate(powi_rf_exp(nx))
  allocate(powe_oh_exp(nx))
  allocate(powe_rad_exp(nx))
  allocate(powe_ion_exp(nx))
  allocate(powi_ion_exp(nx))
  allocate(powe_wdot_exp(nx))
  allocate(powi_wdot_exp(nx))
  allocate(pow_e_fus(nx))
  allocate(pow_i_fus(nx))
  allocate(pow_ei(nx))
  allocate(pow_e(nx))
  allocate(pow_i(nx))
  allocate(powi_cx_exp(nx))
  allocate(flow_wall_exp(nx))
  allocate(flow_beam(nx))
  allocate(flow_mom(nx))
 
  ! Lumped species
  allocate(n_lump_therm(nx))
  allocate(n_lump_fast(nx))

end subroutine allocate_internals
