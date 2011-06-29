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
