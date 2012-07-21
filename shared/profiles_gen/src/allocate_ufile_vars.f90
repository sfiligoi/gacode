subroutine allocate_ufile_vars

  use prgen_read_globals

  implicit none

  allocate(ufile_rho(ufile_nj))
  allocate(ufile_r(ufile_nj))
  allocate(ufile_te(ufile_nj))

end subroutine allocate_ufile_vars
