subroutine allocate_ufile_vars

  use prgen_globals

  implicit none

  allocate(ufile_ne(nx))
  allocate(ufile_nm1(nx))
  allocate(ufile_nm2(nx))
  allocate(ufile_nm3(nx))
  allocate(ufile_te(nx))
  allocate(ufile_ti(nx))
  allocate(ufile_zeff(nx))
  allocate(ufile_pres(nx))
  allocate(ufile_vrot(nx))

end subroutine allocate_ufile_vars
