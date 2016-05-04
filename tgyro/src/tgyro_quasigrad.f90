!------------------------------------------------------------
! tgyro_quasigrad.f90
!
! PURPOSE:
!  Simple routine to enforce derivative of
!
!             ne = sum_i n_i z_i
!
! according to different rules.
!------------------------------------------------------------

subroutine tgyro_quasigrad(ne,dlnnedr,ni,dlnnidr,zi,n_ion)

  use tgyro_globals, only : &
       tgyro_quasineutral_flag, &
       tgyro_fix_concentration_flag, &
       loc_he_feedback_flag, &
       i_ash

  implicit none

  integer, intent(in) :: n_ion
  real, intent(in) :: ne
  real, intent(in) :: dlnnedr
  real, dimension(n_ion), intent(in) :: zi
  real, dimension(n_ion), intent(inout) :: ni
  real, dimension(n_ion), intent(inout) :: dlnnidr

  integer :: i

  if (tgyro_quasineutral_flag == 0) return

  if (n_ion == 1) then

     dlnnidr(1) = dlnnedr

  else

     if (tgyro_fix_concentration_flag == 0) then

        ! Adjust ion 1 gradient to enforce quasineutrality.

        ! Temporary storage 
        dlnnidr(1) = ne*dlnnedr 

        do i=2,n_ion
           dlnnidr(1) = dlnnidr(1)-zi(i)*ni(i)*dlnnidr(i)
        enddo

        dlnnidr(1) = dlnnidr(1)/ni(1)

     else

        ! Adjust all ion gradients at fixed concentration ratios (n2/n1, n3/n1, etc)

        ! Some algebra shows that this implies all gradient scale lengths are equal
        do i=1,n_ion
           if (i /= i_ash .or. loc_he_feedback_flag == 0) dlnnidr(i) = dlnnedr 
        enddo
     
     endif

  endif


end subroutine tgyro_quasigrad

