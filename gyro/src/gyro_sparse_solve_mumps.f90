!---------------------------------
! STUB
!---------------------------------

subroutine gyro_sparse_solve_mumps(n_elem,n_row,matnum,i_solve)

  !---------------------------------------------------
  integer, intent(in) :: n_elem
  integer, intent(in) :: n_row
  integer, intent(in) :: matnum
  integer, intent(in) :: i_solve
  !---------------------------------------------------

  call catch_error('ERROR: (GYRO) MUMPS not available.')

end subroutine gyro_sparse_solve_mumps
