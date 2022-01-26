subroutine allocate_corsica_vars

  use prgen_globals

  implicit none

  allocate(corsica_rho(nx))
  allocate(corsica_r_a(nx))
  allocate(corsica_vl(nx))
  allocate(corsica_ne(nx))
  allocate(corsica_ndt(nx))
  allocate(corsica_nz(nx))
  allocate(corsica_nalpha(nx))
  allocate(corsica_zeff(nx))
  allocate(corsica_j(nx))
  allocate(corsica_jbs(nx))

end subroutine allocate_corsica_vars
