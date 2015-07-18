subroutine cgyro_freq

  use cgyro_globals

  implicit none

  real :: total_weight,dfr,dfi
  real, dimension(n_radial,n_theta) :: mode_weight
  complex, dimension(n_radial,n_theta) :: freq_loc

  if (i_time == 0) return

  if (n == 0) then
     freq = 0.0
     freq_err = 0.0
     return
  endif

  ! Use potential as gauge for frequency
  mode_weight(:,:) = abs(field_old(:,:,1))

  ! Define local frequencies
  if (n > 0) then
     freq_loc(:,:) = (i_c/delta_t)*log(field_old(:,:,1)/field_old2(:,:,1))
  else
     freq_loc = 0.0
  endif

  total_weight = sum(mode_weight(:,:))

  freq = sum(freq_loc(:,:)*mode_weight(:,:))/total_weight

  ! Fractional Frequency Error
  dfr = sum(abs(real(freq_loc(:,:)-freq))*mode_weight(:,:))
  dfi = sum(abs(aimag(freq_loc(:,:)-freq))*mode_weight(:,:))

  freq_err = (dfr+i_c*dfi)/total_weight/abs(freq)

  if (abs(freq_err) < freq_tol) signal=1

end subroutine cgyro_freq
