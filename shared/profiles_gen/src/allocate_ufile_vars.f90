subroutine allocate_ufile_vars

  use prgen_read_globals

  implicit none

  allocate(ufile_ne(nx))
  allocate(ufile_te(nx))
  allocate(ufile_ti(nx))
  allocate(ufile_zeff(nx))

end subroutine allocate_ufile_vars
