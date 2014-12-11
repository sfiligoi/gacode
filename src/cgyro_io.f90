module cgyro_io 

contains

  !-----------------------------------------------------------
  ! cgyro_info.f90
  !
  ! PURPOSE:
  !  Routine to write line to runfile.
  !-----------------------------------------------------------

  subroutine cgyro_info(message)

    use cgyro_globals, only : silent_flag, io, runfile, i_proc, path

    implicit none

    character (len=*), intent(in) :: message

    if (silent_flag == 0 .and. i_proc == 0) then
       open(unit=io,file=trim(path)//runfile,status='old',position='append')
       write(io,'(a)') 'INFO: (CGYRO) '//message
       close(io)
    endif

  end subroutine cgyro_info

  subroutine cgyro_error(message)

    use cgyro_globals, only : error_status, error_message, &
         silent_flag, io, runfile, i_proc, path

    implicit none

    character (len=*), intent(in) :: message

    error_status  = 1
    error_message = message

    if (silent_flag == 0 .and. i_proc == 0) then
       open(unit=io,file=trim(path)//runfile,status='old',position='append')
       write(io,'(a)') 'ERROR: (CGYRO) '//message
       close(io)
    endif

  end subroutine cgyro_error

end module cgyro_io
