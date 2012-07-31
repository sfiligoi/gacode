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
  allocate(ufile_volume(nx))
  allocate(ufile_qnbii(nx))
  allocate(ufile_qnbie(nx))
  allocate(ufile_qicrhi(nx))
  allocate(ufile_qicrhe(nx))
  allocate(ufile_qei(nx))
  allocate(ufile_qrad(nx))
  allocate(ufile_qeche(nx))
  allocate(ufile_qechi(nx))
  allocate(ufile_qohm(nx))
  allocate(ufile_qwalli(nx))
  allocate(ufile_qwalle(nx))

end subroutine allocate_ufile_vars
