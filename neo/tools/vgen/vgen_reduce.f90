subroutine vgen_reduce(vec,n_vec)

  use mpi
  use vgen_globals

  implicit none
  integer, intent(in) :: n_vec
  real, dimension(n_vec), intent(inout) :: vec
  real, dimension(n_vec) :: vec_glob

  call MPI_ALLREDUCE(vec,vec_glob,n_vec, &
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,i_err)

  vec = vec_glob

end subroutine vgen_reduce


