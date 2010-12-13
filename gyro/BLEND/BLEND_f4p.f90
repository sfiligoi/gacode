!-----------------------------------------------------
! BLEND_f4p.f90
!
! PURPOSE:
!  Compute fourth-order blending (cubic b-spline) 
!  function derivative (see BLEND_f4.f90).
!---------------------------------------------

real function BLEND_f4p(m,t0) result (f)

  implicit none

  real, intent(in) :: t0
  integer, intent(in) :: m

  real :: t

  t = t0-m

  if (t <= 0.0 .or. t >= 4.0) then
     f = 0.0 
  else if (t < 1.0) then
     f = 0.75*t*t
  else if (t < 2.0) then
     f = -3.0+6.0*t-2.25*t*t
  else if (t < 3.0) then
     f = 15.0-12.0*t+2.25*t*t
  else
     f = -0.75*(4.0-t)**2
  endif 

  ! Rescale to reflect unit normalization.

  f = (2.0/3.0)*f

end function BLEND_f4p
