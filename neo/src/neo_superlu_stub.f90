subroutine SOLVE_superlu(n_elem, n_size, a, a_iindx, a_jindx)
    implicit none
    integer, intent (in) :: n_elem, n_size
    real, dimension(:), intent (in) :: a
    integer, dimension(:), intent (in) :: a_iindx
    integer, dimension(:), intent (in) :: a_jindx
    call neo_error('ERROR: (NEO) SuperLU package not used')
end subroutine SOLVE_superlu
