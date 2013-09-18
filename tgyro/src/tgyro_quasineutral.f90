!------------------------------------------------------------
! tgyro_quasineutral.f90
!
! PURPOSE:
!  Simple routine to enforce quasineutrality.  This is only
!  used by the tgyro_global_* mode.
!------------------------------------------------------------

subroutine tgyro_quasineutral(ni,ne,dlnnidr,dlnnedr,zi,n_ion,n_r)

  use tgyro_globals, only : tgyro_quasineutral_flag

  implicit none

  integer, intent(in) :: n_ion
  integer, intent(in) :: n_r
  real, intent(inout) :: ni(n_ion,n_r)
  real, intent(in) :: ne(n_r)
  real, intent(inout) :: dlnnidr(n_ion,n_r)
  real, intent(in) :: dlnnedr(n_r)
  real, intent(in) :: zi(n_ion)

  integer :: i

  if (tgyro_quasineutral_flag == 0) return
 
  if (n_ion == 1) then

     ! ni = ne
     ni(1,:) = ne(:)
     dlnnidr(1,:) = dlnnedr(:)

  else

     ! ni = ne-z*nz 

     ni(1,:) = ne(:)
     dlnnidr(1,:) = ne(:)*dlnnedr(:) 

     do i=2,n_ion

        ni(1,:) = ni(1,:)-zi(i)*ni(i,:)
        dlnnidr(1,:) = dlnnidr(1,:)-zi(i)*ni(i,:)*dlnnidr(i,:)

     enddo

     dlnnidr(1,:) = dlnnidr(1,:)/ni(1,:)

  endif

end subroutine tgyro_quasineutral
