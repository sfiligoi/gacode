subroutine allocate_peqdsk_vars

  use prgen_globals

  implicit none

  allocate(peqdsk_psi(peqdsk_nj))
  allocate(peqdsk_ne(peqdsk_nj))
  allocate(peqdsk_te(peqdsk_nj))
  allocate(peqdsk_ni(peqdsk_nj))
  allocate(peqdsk_ti(peqdsk_nj))
  allocate(peqdsk_nb(peqdsk_nj))
  allocate(peqdsk_pb(peqdsk_nj))
  allocate(peqdsk_omegat(peqdsk_nj))
  allocate(peqdsk_omgeb(peqdsk_nj))

end subroutine allocate_peqdsk_vars
