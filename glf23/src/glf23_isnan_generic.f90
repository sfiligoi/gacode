  logical function glf23_isnan(x)
  ! 
  ! test function for platforms where isnan is not intrinsic.
  ! This numerical test may not work depending on the compiler
  !
  implicit none
  real, intent(in) :: x

  glf23_isnan = .false.

  if (x/=x) glf23_isnan = .true.

  end function glf23_isnan
   
      
