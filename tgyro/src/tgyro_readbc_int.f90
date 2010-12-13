subroutine tgyro_readbc_int(p)

  use tgyro_globals

  implicit none
  
  integer, intent(inout) :: p
  integer :: i_err

  include 'mpif.h'

  if (i_proc_global == 0) read(1,*) p

  call MPI_BCAST(p,&
       1 ,&
       MPI_INTEGER,&
       0 ,&
       MPI_COMM_WORLD,&
       i_err)

end subroutine tgyro_readbc_int
