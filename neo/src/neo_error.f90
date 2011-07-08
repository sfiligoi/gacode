subroutine neo_error(message)

  use neo_globals, only : error_status, error_message, write_out_mode

  implicit none
  character (len=*), intent(in) :: message

  error_status  = 1
  error_message = message

  if (write_out_mode > 1) then
     print '(a)', message
  endif

end subroutine neo_error
