subroutine tglf_error(itype,error_message)
  
  ! error handeling routine for TGLF
  ! simply writes the error or informatin message
  ! if itype = 1 the error was fatal so TGLF is shutdown
  
  implicit none
  integer, intent(in):: itype
  character(len=32), intent(in) :: error_message
 
  print *,'ERROR: (tglf) ',error_message

  if (itype == 1) call tglf_shutdown
 
end subroutine tglf_error

logical function tglf_isnan(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  tglf_isnan = .false.

  if (x /= x) tglf_isnan = .true.

end function tglf_isnan

logical function tglf_isinf(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  tglf_isinf = .false.

  if (0.0*x /= 0.0) tglf_isinf = .true.

end function tglf_isinf
   
      
