subroutine neo_scalapack(info)

  implicit none

  integer, intent(inout) :: info

  call neo_error('ERROR: (NEO) SCALAPACK library not available.')

  info = 1

end subroutine neo_scalapack
