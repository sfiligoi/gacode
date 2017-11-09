module neo_sparse_solve

  implicit none

  public :: SOLVE_sparse

contains

  subroutine SOLVE_sparse(n_elem, n_size, a, a_iindx, a_jindx)
    use neo_globals
    use neo_umfpack
    use neo_superlu
    implicit none
    integer, intent (in) :: n_elem, n_size
    real, dimension(:), intent (in) :: a
    integer, dimension(:), intent (in) :: a_iindx
    integer, dimension(:), intent (in) :: a_jindx

#ifdef USE_UMFPACK
    call SOLVE_umfpack(n_elem, n_size, a, a_iindx, a_jindx)
#else
    call SOLVE_superlu(n_elem, n_size, a, a_iindx, a_jindx)
#endif
  end subroutine SOLVE_sparse

end module neo_sparse_solve
