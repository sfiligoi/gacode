!-----------------------------------------------------
! BLEND_F.f90
!
! PURPOSE:
!  Evaluate blending function basis vectors.  
!------------------------------------------------------

complex function BLEND_F(m,x,p_IN)

  use BLEND_private

  implicit none

  integer, intent(in) :: m
  real, intent(in) :: x
  complex, intent(in) :: p_IN

  real :: t

  real, external :: BLEND_f2
  real, external :: BLEND_f3
  real, external :: BLEND_f4


  t = (x+1.0)/d

  ! EVALUATION OF computed blending function:
  !
  ! use f = fn(m)+phase*fn(m+n_fit) 

  select case (n)

  case (2)

     BLEND_F = BLEND_f2(m-n,t)+p_IN*BLEND_f2(m-n+n_fit,t)

  case (3)

     BLEND_F = BLEND_f3(m-n,t)+p_IN*BLEND_f3(m-n+n_fit,t)

  case (4)

     BLEND_F = BLEND_f4(m-n,t)+p_IN*BLEND_f4(m-n+n_fit,t)

  end select

end function BLEND_F
