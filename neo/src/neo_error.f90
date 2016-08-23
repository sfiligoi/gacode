subroutine neo_error(message)

  use neo_globals, only : error_status, error_message, &
       silent_flag, io_neoout, runfile_neoout, i_proc, path

  implicit none
  character (len=*), intent(in) :: message

  error_status  = 1
  error_message = message

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
        write(io_neoout,'(t2,a)') message
        close(io_neoout)
  endif

end subroutine neo_error
