!-----------------------------------------------------
! get_ampere_explicit
!
! PURPOSE:
!  Solve explicit field equation for Apar.
!---------------------------------------------

subroutine get_ampere_explicit

  use gyro_globals
  use gyro_collision_private

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  ! Generate RHS (vel_sum_a) for Apar part of field-solve.
  call gyro_velocity_sum(2)

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_ampere,n_ampere_row,2,1)
  else
     call sparse_solve_mumps(n_ampere,n_ampere_row,2,1)
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[get_ampere_explicit done]'
  endif

end subroutine get_ampere_explicit
