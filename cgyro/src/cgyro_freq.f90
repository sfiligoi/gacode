!---------------------------------------------------------
! cgyro_freq.f90
!
! PURPOSE:
!  Compute estimates of linear eigenvalue w where:
!
!                   -iwt
!            phi ~ e
!---------------------------------------------------------

subroutine cgyro_freq

  use cgyro_globals

  implicit none

  real :: total_weight,dfr,dfi
  real, dimension(nc) :: mode_weight
  complex, dimension(nc,my_toroidal:my_toroidal) :: freq_loc

  if (i_time == 0) then

     freq = 0.0
     freq_err = 0.0

  else

     !--------------------------------------------------
     ! Standard method: sum all wavenumbers at a given n
     !--------------------------------------------------

     ! Use potential to compute frequency
     mode_weight(:) = abs(field_old(1,:,my_toroidal))

     ! Define local frequencies
     do ic=1,nc
        if (abs(field_old(1,ic,my_toroidal)) > 1e-12 .and. abs(field_old2(1,ic,my_toroidal)) > 1e-12) then
           freq_loc(ic,my_toroidal) = (i_c/delta_t)*log(field_old(1,ic,my_toroidal)/field_old2(1,ic,my_toroidal))
        else
           freq_loc(ic,my_toroidal) = 0.0
        endif
     enddo

     total_weight = sum(mode_weight(:))

     freq = sum(freq_loc(:,my_toroidal)*mode_weight(:))/total_weight

     ! Fractional Frequency Error
     dfr = sum(abs(real(freq_loc(:,my_toroidal)-freq))*mode_weight(:))
     dfi = sum(abs(aimag(freq_loc(:,my_toroidal)-freq))*mode_weight(:))

     freq_err = (dfr+i_c*dfi)/total_weight/abs(freq)

     if (n_toroidal == 1 .and. abs(freq_err) < freq_tol) signal=1

  endif


end subroutine cgyro_freq
