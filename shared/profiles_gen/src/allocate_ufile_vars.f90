subroutine allocate_ufile_vars

  use prgen_read_globals

  implicit none

  allocate(ufile_rho(ufile_nj))
  allocate(ufile_te(ufile_nj))
  allocate(ufile_ti(ufile_nj))
  allocate(ufile_zeff(ufile_nj))
  allocate(ufile_ene(ufile_nj))
  allocate(ufile_en(ufile_nj,ufile_nion))
  allocate(ufile_enbeam(ufile_nj))
  allocate(ufile_pfast(ufile_nj))
  allocate(ufile_ptot(ufile_nj))
  allocate(ufile_angrot(ufile_nj))

end subroutine allocate_ufile_vars
