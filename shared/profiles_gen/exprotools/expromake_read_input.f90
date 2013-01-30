subroutine expromake_read_input

  use expromake_globals

  implicit none

  open(unit=1,file='input.expromake.gen',status='old')
  read(1,*) z(1)
  read(1,*) z(2)
  read(1,*) z(3)
  close(1)

end subroutine expromake_read_input
