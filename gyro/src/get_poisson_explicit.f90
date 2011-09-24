!-----------------------------------------------------
! get_poisson_explicit
!
! PURPOSE:
!  Solve explicit field equation for phi.
!
! REVISIONS
! 08 Sept 06: jc
!  Created.
!---------------------------------------------

subroutine get_poisson_explicit

  use gyro_globals
  use gyro_poisson_private

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  ! Generate RHS for phi part of field-solve.
  call gyro_velocity_sum(1)

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_poisson,n_poisson_row,1,1)
  else
     call sparse_solve_mumps(n_poisson,n_poisson_row,1,1)
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[get_poisson_explicit done]'
  endif

end subroutine get_poisson_explicit
