subroutine gacode_system(arg)

  implicit none
  character(len=*), intent(in) :: arg

  call execute_command_line(arg)

end subroutine gacode_system
