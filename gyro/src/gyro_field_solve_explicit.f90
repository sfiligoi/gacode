!-----------------------------------------------------
! gyro_field_solve_explicit.f90
!
! PURPOSE:
!  This routine coordinates calling of the appropriate 
!  sequence of routines for solution of the EXPLICIT 
!  Maxwell equations.
!-----------------------------------------------------

subroutine gyro_field_solve_explicit

  use gyro_globals

  implicit none

  call gyro_operators_on_h

  ! Velocity-space sums (RHS)
  do ix=1,n_field
     call gyro_velocity_sum(ix)
  enddo

  call gyro_timer_in('Field-explicit')

  select case (n_field)

  case (1)

     if (sparse_method == 1) then
        call sparse_solve_umfpack(n_poisson,n_poisson_row,1,1)
     else
        call sparse_solve_mumps(n_poisson,n_poisson_row,1,1)
     endif

  case (2)

     if (sparse_method == 1) then
        call sparse_solve_umfpack(n_poisson,n_poisson_row,1,1)
        call sparse_solve_umfpack(n_ampere,n_ampere_row,2,1)
     else
        call sparse_solve_mumps(n_poisson,n_poisson_row,1,1)
        call sparse_solve_mumps(n_ampere,n_ampere_row,2,1)
     endif

  case (3)

     if (sparse_method == 1) then
        call sparse_solve_umfpack(n_poissonaperp,n_poissonaperp_row,4,1)
        call sparse_solve_umfpack(n_ampere,n_ampere_row,2,1)
     else
        call sparse_solve_mumps(n_poissonaperp,n_poissonaperp_row,4,1)
        call sparse_solve_mumps(n_ampere,n_ampere_row,2,1)
     endif

  end select

  call gyro_timer_out('Field-explicit')

  call gyro_field_interpolation

end subroutine gyro_field_solve_explicit
