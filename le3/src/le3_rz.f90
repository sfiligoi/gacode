subroutine le3_rz(t0,p0,r0,z0,drdt0,drdp0,dzdt0,dzdp0,jac0,rc0)

  use le3_globals

  implicit none

  real, intent(in) :: t0,p0
  real, intent(inout) :: r0,z0,jac0
  real, intent(inout) :: drdt0,dzdt0
  real, intent(inout) :: drdp0,dzdp0
  real, intent(inout) :: rc0
  real :: drdr0,dzdr0
  real :: drdtt0,dzdtt0
  real :: ang,x,a,a_t,a_tt

  ang = m*t0+n*p0

  x    = asin(delta)
  a    = t0+x*sin(t0)
  a_t  = 1.0+x*cos(t0)
  a_tt = -x*sin(t0)

  ! R
  r0 = rmaj + rmin*cos(a) + hmin*cos(ang)

  ! dR/dr
  drdr0 = shift+cos(a)-s_delta/cos(x)*sin(t0)*sin(a) + dhmindr*cos(ang)

  ! dR/d(tb)
  drdt0 = -rmin*sin(a)*a_t - m*hmin*sin(ang) 

  ! dR/d(pb)
  drdp0 = -n*hmin*sin(ang)                 

  ! d^2R/d(tb)^2
  drdtt0 = -rmin*a_t**2*cos(a)-rmin*a_tt*sin(a)-m**2*hmin*cos(ang)

  a    = t0+zeta*sin(2*t0)
  a_t  = 1.0+2.0*zeta*cos(2*t0)
  a_tt = -4*zeta*sin(2*t0)

  ! Z
  z0 = zmag+kappa*rmin*sin(a) + hmin*sin(ang)

  ! dZ/dr
  dzdr0 = dzmag+kappa*(1.0+s_kappa)*sin(a)+kappa*s_zeta*cos(a)*sin(2*t0) &
       + dhmindr*sin(ang)

  ! dZ/d(tb)
  dzdt0 =  kappa*rmin*cos(a)*a_t + m*hmin*cos(ang) 

  ! dZ/d(pb)
  dzdp0 = n*hmin*cos(ang)

  ! d^2Z/d(tb)^2
  dzdtt0 = -kappa*rmin*sin(a)*a_t**2+kappa*rmin*cos(a)*a_tt

  ! J [dtb/dr terms vanish]
  jac0 = r0*(drdr0*dzdt0-drdt0*dzdr0)

  ! radius of curvature r_c (this is independent of theta coordinates).
  rc0 = (drdt0**2+dzdt0**2)**1.5/(drdt0*dzdtt0-dzdt0*drdtt0)  
 
end subroutine le3_rz
