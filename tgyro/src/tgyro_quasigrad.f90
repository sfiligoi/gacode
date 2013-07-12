!------------------------------------------------------------
! tgyro_quasigrad.f90
!
! PURPOSE:
!  Simple routine to enforce derivative of
!
!             ne = sum_i n_i z_i
!
! by evolving 
!------------------------------------------------------------

subroutine tgyro_quasigrad(ne,dlnnedr,ni,dlnnidr,zi,n_ion)

  use tgyro_globals, only : tgyro_fix_concentration_flag 

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

     if (tgyro_fix_concentration_flag == 0) then

        ! Temporary storage 
        dlnnidr(1) = ne*dlnnedr 

        do i=2,n_ion
           dlnnidr(1) = dlnnidr(1)-zi(i)*ni(i)*dlnnidr(i)
        enddo

        dlnnidr(1) = dlnnidr(1)/ni(1)

     endif

  endif


end subroutine tgyro_quasigrad

