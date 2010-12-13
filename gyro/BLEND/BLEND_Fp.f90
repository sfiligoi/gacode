!-----------------------------------------------------
! BLEND_Fp.f90
!
! PURPOSE:
!  Evaluate blending function derivative.  
!------------------------------------------------------

complex function BLEND_Fp(m,x,p_IN) result (f)

  use BLEND_private

  implicit none

  complex, intent(in) :: p_IN

  real, intent(in) :: x
  integer, intent(in) :: m

  real :: t
  real :: a

  real, external :: BLEND_f2p
  real, external :: BLEND_f3p
  real, external :: BLEND_f4p


  t = (x+1.0)/d

  ! EVALUATION OF computed blending function:
  !
  ! use f = fn(m)+phase*fn(m+n_fit) 

  select case (n)

  case (2)

     f = BLEND_f2p(m-n,t)+p_IN*BLEND_f2p(m-n+n_fit,t)

  case (3)

     f = BLEND_f3p(m-n,t)+p_IN*BLEND_f3p(m-n+n_fit,t)

  case (4)

     f = BLEND_f4p(m-n,t)+p_IN*BLEND_f4p(m-n+n_fit,t)

  case default

     ! This code is meant to cause GYRO to fail.
     ! This is done, rather than simply calling 
     ! catch_error(), as a workaround to enable
     ! the X1 to vectorize.

     a = 0.0
     a = 1.0/a
     t = 1.0/a
     f = a*t

  end select

  ! Scale derivative (d/dx)

  f = f/d

end function BLEND_Fp
