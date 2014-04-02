subroutine tglf_error(itype,error_message)

  
  ! error handeling routine for TGLF
  ! simply writes the error or informatin message
  ! if itype = 1 the error was fatal so TGLF is shutdown
  
  implicit none
  integer, intent(in):: itype
  character(len=*), intent(in) :: error_message
 
  print *,'ERROR: (tglf) ',error_message

  if (itype == 1) call tglf_shutdown
 
end subroutine tglf_error

logical function tglf_isnan(x)
  ! 
  ! test real variables for NAN 
  !
  use ieee_arithmetic, only: ieee_is_nan
  implicit none
  real, intent(in) :: x

  tglf_isnan  = ieee_is_nan(x)

end function tglf_isnan

logical function tglf_isinf(x)
  ! 
  ! test real variables for Inf
  !
  implicit none
  real, intent(in) :: x

  tglf_isinf = (ABS(x) > HUGE(x))

end function tglf_isinf
   
      
