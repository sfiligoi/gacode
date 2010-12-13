subroutine tgyro_quasigrad(ne,dlnnedr,ni,dlnnidr,zi,n_ion)

  implicit none

  integer, intent(in) :: n_ion
  real, intent(in) :: ne
  real, intent(in) :: dlnnedr
  real, dimension(n_ion), intent(in) :: zi
  real, dimension(n_ion), intent(inout) :: ni
  real, dimension(n_ion), intent(inout) :: dlnnidr

  integer :: i

  if (n_ion == 1) then

     dlnnidr(1) = dlnnedr

  else

     dlnnidr(1) = ne*dlnnedr 

     do i=2,n_ion
        dlnnidr(1) = dlnnidr(1)-zi(i)*ni(i)*dlnnidr(i)
     enddo

     dlnnidr(1) = dlnnidr(1)/ni(1)

  endif


end subroutine tgyro_quasigrad

