subroutine cgyro_freq

  use cgyro_globals

  implicit none

  real :: total_weight,dfr,dfi
  real, dimension(n_radial,n_theta) :: mode_weight
  complex, dimension(n_radial,n_theta) :: freq_loc

  if (itime == 0) return

  ! Use potential as gauge for frequency
  mode_weight(:,:) = abs(field(:,:,1))

  ! Define local frequencies
  freq_loc(:,:) = (i_c/delta_t)*log(field(:,:,1)/field_old(:,:,1))

  total_weight = sum(mode_weight(:,:))

  freq = sum(freq_loc(:,:)*mode_weight(:,:))/total_weight

  ! Fractional Frequency Error
  dfr = sum(abs(real(freq_loc(:,:)-freq))*mode_weight(:,:))
  dfi = sum(abs(imag(freq_loc(:,:)-freq))*mode_weight(:,:))

  freq_err = (dfr+i_c*dfi)/total_weight/abs(freq)

  if (abs(freq_err) < freq_tol) signal=1

end subroutine cgyro_freq
