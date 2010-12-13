!-----------------------------------------------------
! get_maxwell3_explicit
!
! PURPOSE:
!  Solve explicit field equation for phi, apar, and aperp
!
!---------------------------------------------

subroutine get_maxwell3_explicit

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  ! Get A_parallel solution (uncoupled field solve)
  call get_ampere_explicit

  ! Get Phi, B_parallel solutions (coupled field equations)

  ! Generate RHS (vel_sum_p) for phi part of field-solve.
  ! uses gyro_h = (G h)=J0 h
  call get_vel_sum_p
  
  ! Generate RHS (vel_sum_aperp) for bpar part of field-solve.
  ! uses gyro_h_perp = (G_aperp h)=1/2(J0+J1) h
  call get_gyro_h_aperp
  call get_vel_sum_aperp

  if (sparse_method == 1) then
     call sparse_solve_umfpack(n_poissonaperp,n_poissonaperp_row,4,1)
  else
     call sparse_solve_mumps(n_poissonaperp,n_poissonaperp_row,4,1)
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[get_maxwell3_explicit done]'
  endif

end subroutine get_maxwell3_explicit
