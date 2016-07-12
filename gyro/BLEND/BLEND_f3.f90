!-----------------------------------------------------
! BLEND_f3.f90
!
! PURPOSE:
!  Compute third-order blending function.
! 
! NOTES:
!                 2
!  f(m,t0) = 0.5 t                 if 0 < t < 1
!
!                     2
!          = 0.5 (-2 t  + 6t - 3)  if 1 < t < 2
!
!                     2
!          = 0.5 (3-t)             if 2 < t < 3
!
!  where t = t0-m
!---------------------------------------------

real function BLEND_f3(m,t0)

  implicit none

  real, intent(in) :: t0
  integer, intent(in) :: m

  real :: t

  t = t0-m

  if (t <= 0.0 .or. t >= 3.0) then
     BLEND_f3 = 0.0 
  else if (t < 1.0) then
     BLEND_f3 = 0.5*t*t
  else if (t < 2.0) then
     BLEND_f3 = 0.5*((-2.0*t+6.0)*t-3.0)
  else
     BLEND_f3 = 0.5*(3.0-t)**2
  endif

end function BLEND_f3
