  logical function gftm_isnan(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  gftm_isnan = ISNAN(x)

  end function gftm_isnan
   
      
