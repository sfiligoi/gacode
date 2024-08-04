  SUBROUTINE gftm_error(itype,error_message)
  
  ! error handeling routine for gftm
  ! simply writes the error or informatin message
  ! if itype = 1 the error was fatal so gftm is shutdown
  
  implicit none
  integer, intent(in):: itype
  character(len=*), intent(in) :: error_message
 
  print *,'ERROR: (gftm) ',error_message

  if (itype == 1) then
     CALL gftm_shutdown
     STOP
  endif
 
  END SUBROUTINE gftm_error   
      
