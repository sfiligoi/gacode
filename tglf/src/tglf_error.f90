  SUBROUTINE tglf_error(itype,error_message)
  
  ! error handeling routine for TGLF
  ! simply writes the error or informatin message
  ! if itype = 1 the error was fatal so TGLF is shutdown
  
  implicit none
  integer, intent(in):: itype
  character(len=*), intent(in) :: error_message
 
  print *,'ERROR: (tglf) ',error_message

  if (itype == 1) then
     CALL tglf_shutdown
     STOP
  endif
 
  END SUBROUTINE tglf_error

  logical function tglf_isnan(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  tglf_isnan = .false.

  if (ISNAN(x)) tglf_isnan = .true.

  end function tglf_isnan

  logical function tglf_isinf(x)
  ! 
  ! test real variables for NAN 
  !
  implicit none
  real, intent(in) :: x

  tglf_isinf = .false.

  if (ABS(x) > HUGE(x)) tglf_isinf = .true.

  end function tglf_isinf
   
      
