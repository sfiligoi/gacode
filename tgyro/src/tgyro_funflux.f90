real function tgyro_funflux(r,zi,ze,ix,is)

  implicit none

  real, intent(in) :: r,ix,is
  real, intent(in) :: zi,ze

  tgyro_funflux = 0.2*((zi-1.130)/(1.136-1.130)+(ze-1.226)/(1.471-1.226))

end function tgyro_funflux
