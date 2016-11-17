!-----------------------------------------------------
! BLEND_f2.f90
!
! PURPOSE:
!  Compute second-order blending (tent) function.
! 
! NOTES:
!  
!  f(m,t0) = t     if 0 < t < 1
!          = 2-t   if 1 < t < 2
!
!  where t = t0-m.
!---------------------------------------------

real function BLEND_f2(m,t0)

  implicit none

  real, intent(in) :: t0
  integer, intent(in) :: m

  real :: t

  t = t0-m

  if (t <= 0.0 .or. t >= 2.0) then
     BLEND_f2 = 0.0 
  else if (t < 1.0) then
     BLEND_f2 = t
  else 
     BLEND_f2 = 2.0-t
  endif

end function BLEND_f2
