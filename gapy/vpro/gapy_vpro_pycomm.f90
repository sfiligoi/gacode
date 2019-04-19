subroutine vpro_icomm(p)

  implicit none

  integer, intent(inout) :: p

  read(1,*) p

end subroutine vpro_icomm

subroutine vpro_rcomm(x)

  implicit none

  double precision, intent(inout) :: x

  read(1,10) x

10 format(1pe14.7)

end subroutine vpro_rcomm

subroutine vpro_acomm(x,n)

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x

  read(1,10) x

10 format(1pe14.7)

end subroutine vpro_acomm
