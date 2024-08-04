  logical function gftm_isnan(x)
  ! 
  ! test real variables for NAN 
  !
  use ieee_arithmetic, only : ieee_is_nan
  implicit none
  real, intent(in) :: x

  gftm_isnan = IEEE_IS_NAN(x)

  end function gftm_isnan
   
      
