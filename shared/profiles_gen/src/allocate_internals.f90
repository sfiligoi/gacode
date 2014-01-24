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

  allocate(pow_e_nb(nx))
  allocate(pow_i_nb(nx))
  allocate(pow_e_rf(nx))
  allocate(pow_i_rf(nx))
  allocate(pow_e_ohm(nx))
  allocate(powe_ion_exp(nx))
  allocate(powi_ion_exp(nx))
  allocate(powe_wdot_exp(nx))
  allocate(powi_wdot_exp(nx))
  allocate(pow_e_fus(nx))
  allocate(pow_i_fus(nx))
  allocate(pow_ei(nx))
  allocate(pow_e(nx))
  allocate(pow_i(nx))
  allocate(pow_e_rad(nx))
  allocate(pow_e_sync(nx))
  allocate(pow_e_brem(nx))
  allocate(pow_e_line(nx))
  allocate(pow_e_aux(nx))
  allocate(pow_i_aux(nx))
  allocate(powi_cx_exp(nx))
  allocate(flow_wall_exp(nx))
  allocate(flow_beam(nx))
  allocate(flow_mom(nx))

  ion_lookup(:) = '?'
  ion_lookup(1) = 'H'
  ion_lookup(2) = 'D'
  ion_lookup(3) = 'T'
  ion_lookup(4) = 'He4'
  ion_lookup(9) = 'Be'
  ion_lookup(12) = 'C'
  ion_lookup(40) = 'Ar'

end subroutine allocate_internals
