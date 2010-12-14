!-----------------------------------------------------
! get_poisson_solution
!
! PURPOSE:
!  Solve implicit field equation for phi.
!---------------------------------------------

subroutine get_poisson_solution

  use gyro_globals
  use gyro_maxwell_private

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  ! Generate RHS (vel_sum_p) for phi part of field-solve.
  call get_vel_sum_p

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_maxwell,n_maxwell_row,3,1)
  else
     call sparse_solve_mumps(n_maxwell,n_maxwell_row,3,1)
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[get_poisson_solution done]'
  endif

end subroutine get_poisson_solution
