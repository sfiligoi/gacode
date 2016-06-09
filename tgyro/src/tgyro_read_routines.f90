subroutine tgyro_readbc_real(xr)

  use mpi
  use tgyro_globals

  implicit none

  real, intent(inout) :: xr
  integer :: i_err


  if (i_proc_global == 0) read(1,*) xr

  call MPI_BCAST(xr,&
       1 ,&
       MPI_DOUBLE_PRECISION,&
       0 ,&
       MPI_COMM_WORLD,&
       i_err)

end subroutine tgyro_readbc_real

subroutine tgyro_readbc_int(p)

  use mpi
  use tgyro_globals

  implicit none
  
  integer, intent(inout) :: p
  integer :: i_err

  if (i_proc_global == 0) read(1,*) p

  call MPI_BCAST(p,&
       1,&
       MPI_INTEGER,&
       0,&
       MPI_COMM_WORLD,&
       i_err)

end subroutine tgyro_readbc_int

subroutine tgyro_readbc_char(xc)

  use mpi
  use tgyro_globals

  implicit none
  
  character(len=5), intent(inout) :: xc
  integer :: i_err

  if (i_proc_global == 0) read(1,*) xc

  call MPI_BCAST(xc,&
       5,&
       MPI_CHARACTER,&
       0,&
       MPI_COMM_WORLD,&
       i_err)

end subroutine tgyro_readbc_char
