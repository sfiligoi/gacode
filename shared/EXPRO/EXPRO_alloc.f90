subroutine EXPRO_alloc(path,flag)

  use EXPRO_interface

  implicit none

  character(len=*), intent(in) :: path
  integer, intent(in) :: flag

  call EXPRO_alloc_control(0,path,flag)

end subroutine EXPRO_alloc
