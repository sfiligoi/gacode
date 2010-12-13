!-----------------------------------------------------
! BLEND_f4.f90
!
! PURPOSE:
!  Compute fourth-order blending (cubic b-spline) 
!  function 
!  
! NOTES:
!                        3
!  (3/2) f(m,t0) = 0.25 t                 if 0 < t < 1
!
!                        2         3
!          = 1 - 3t + 3 t  - 0.75 t       if 1 < t < 2
!
!                            2         3
!          = -11 + 15 t - 6 t  + 0.75 t   if 2 < t < 3
!
!                      3
!          = 0.25 (4-t)                   if 3 < t < 4
! 
!  where t = t0-m
!---------------------------------------------

real function BLEND_f4(m,t0) result (f)

  implicit none

  real, intent(in) :: t0
  integer, intent(in) :: m

  real :: t

  t = t0-m

  if (t <= 0.0 .or. t >= 4.0) then
     f = 0.0 
  else if (t < 1.0) then
     f = t*t*t/6.0
  else if (t < 2.0) then
     f = (1.0+t*(-3.0+t*(3.0-0.75*t)))*2.0/3.0
  else if (t < 3.0) then
     f = (-11.0+t*(15.0+t*(-6.0+0.75*t)))*2.0/3.0
  else
     f = (4.0-t)**3/6.0
  endif

end function BLEND_f4
