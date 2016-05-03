subroutine gacode_system(arg)

  implicit none
  character(len=*), intent(in) :: arg

  call system(arg)

end subroutine gacode_system
