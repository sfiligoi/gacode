!-----------------------------------------------------
! BLEND_f3p.f90
!
! PURPOSE:
!  Compute third-order blending function derivatives
!  (see BLEND_f3.f90).
!-----------------------------------------------------

real function BLEND_f3p(m,t0) result (f)

  implicit none

  real, intent(in) :: t0
  integer, intent(in) :: m

  real :: t

  t = t0-m

  if (t <= 0.0 .or. t >= 3.0) then
     f = 0.0 
  else if (t < 1.0) then
     f = t
  else if (t < 2.0) then
     f = 0.5*(-4.0*t+6.0)
  else
     f = t-3.0
  endif

end function BLEND_f3p
