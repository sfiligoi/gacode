subroutine gkcoll_error(message)

  use gkcoll_globals, only : error_status, error_message, &
       silent_flag, io_gkcollout, runfile_gkcollout, i_proc, path

  implicit none
  character (len=*), intent(in) :: message

  error_status  = 1
  error_message = message

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_gkcollout,file=trim(path)//runfile_gkcollout,&
          status='old',position='append')
        write(io_gkcollout,'(a)') message
        close(io_gkcollout)
  endif

end subroutine gkcoll_error
