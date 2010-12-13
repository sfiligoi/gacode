!-----------------------------------------------------
! get_maxwell_solution
!
! PURPOSE:
!  Solve implicit field equation for phi+Apar.
!---------------------------------------------

subroutine get_maxwell_solution

  use gyro_globals
  use maxwell_private

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  ! Generate RHS (vel_sum_p) for phi part of field-solve.
  call get_vel_sum_p

  ! Generate RHS (vel_sum_a) for Apar part of field-solve.
  call get_vel_sum_a

  if(n_field == 3) then
     call get_gyro_h_aperp
     call get_vel_sum_aperp
  endif

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_maxwell,n_maxwell_row,3,1)
  else
     call sparse_solve_mumps(n_maxwell,n_maxwell_row,3,1)
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[get_maxwell_solution done]'
  endif

end subroutine get_maxwell_solution
