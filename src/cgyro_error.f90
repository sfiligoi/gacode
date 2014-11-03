subroutine cgyro_error(message)

  use cgyro_globals, only : error_status, error_message, &
       silent_flag, io_cgyroout, runfile, i_proc, path

  implicit none
  character (len=*), intent(in) :: message

  error_status  = 1
  error_message = message

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_cgyroout,file=trim(path)//runfile,&
          status='old',position='append')
        write(io_cgyroout,'(a)') message
        close(io_cgyroout)
  endif

end subroutine cgyro_error
