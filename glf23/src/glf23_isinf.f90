
  logical function glf23_isinf(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  glf23_isinf = .false.

  if (ABS(x) > HUGE(x)) glf23_isinf = .true.

  end function glf23_isinf
   
      
