
  logical function gftm_isinf(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  gftm_isinf = .false.

  if (ABS(x) > HUGE(x)) gftm_isinf = .true.

  end function gftm_isinf
   
      
