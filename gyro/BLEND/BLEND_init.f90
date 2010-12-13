!-----------------------------------------------------
! BLEND_init.f90
!
! PURPOSE:
!  Initialize BLEND machinery.
!---------------------------------------------

subroutine BLEND_init(n_IN,n_fit_IN)

  use BLEND_private

  implicit none

  integer, intent(in) :: n_IN
  integer, intent(in) :: n_fit_IN
 
  ! Copy n_IN (order of blend function) and n_fit_IN 
  ! (number of blend functions) to internal variable n.

  n     = n_IN
  n_fit = n_fit_IN

  if (n_fit < 2*n-1) then
     print *,'BLEND_init: not enough fit points.'
     stop
  endif

  d = 2.0/n_fit

  allocate(i_piv(n_fit))

  ! Make blending matrix

  allocate(cs(n_fit,n_fit))
  allocate(c0(n_fit))

end subroutine BLEND_init
