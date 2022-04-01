  logical function gftm_isnan(x)
  ! 
  ! test function for platforms where isnan is not intrinsic.
  ! This numerical test may not work depending on the compiler
  !
  implicit none
  real, intent(in) :: x

  gftm_isnan = .false.

  if (x/=x) gftm_isnan = .true.

  end function gftm_isnan
   
      
