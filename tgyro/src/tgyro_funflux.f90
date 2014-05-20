real function tgyro_funflux(r,zi,ze,ix,is)

  implicit none

  real, intent(in) :: r,ix,is
  real, intent(in) :: zi,ze

  tgyro_funflux = 0.1+zi**2+ze**2

end function tgyro_funflux
