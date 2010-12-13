!-----------------------------------------------------
! BLEND_f2p.f90
!
! PURPOSE:
!  Compute derivative of second-order blending 
!  (tent) function (see BLEND_f2.f90).
!---------------------------------------------

real function BLEND_f2p(m,t0) result (f)

  implicit none

  real, intent(in) :: t0
  integer, intent(in) :: m

  real :: t

  t = t0-m

  if (t <= 0.0 .or. t >= 2.0) then
     f = 0.0 
  else if (t < 1.0) then
     f = 1.0
  else 
     f = -1.0
  endif

end function BLEND_f2p
