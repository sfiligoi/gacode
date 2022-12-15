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
  integer :: itor
  complex, dimension(nc) :: freq_loc

  if (i_time == 0) then

    freq(:) = 0.0
    freq_err(:) = 0.0

  else

    !--------------------------------------------------
    ! Standard method: sum all wavenumbers at a given n
    !--------------------------------------------------

    do itor=nt1,nt2
     ! Use potential to compute frequency
     ! NOTE: Do it once per itor
     mode_weight(:) = abs(field_old(1,:,itor))

     ! Define local frequencies
     do ic=1,nc
        if (abs(field_old(1,ic,itor)) > 1e-12 .and. abs(field_old2(1,ic,itor)) > 1e-12) then
           freq_loc(ic) = (i_c/delta_t)*log(field_old(1,ic,itor)/field_old2(1,ic,itor))
        else
           freq_loc(ic) = 0.0
        endif
     enddo

     total_weight = sum(mode_weight(:))
     freq(itor) = sum(freq_loc(:)*mode_weight(:))/total_weight

     ! Fractional Frequency Error
     dfr = sum(abs(real(freq_loc(:)-freq(itor)))*mode_weight(:))
     dfi = sum(abs(aimag(freq_loc(:)-freq(itor)))*mode_weight(:))
     freq_err(itor) = (dfr+i_c*dfi)/total_weight/abs(freq(itor))

     if (n_toroidal == 1 .and. abs(freq_err(itor)) < freq_tol) signal=1
    enddo

  endif


end subroutine cgyro_freq
