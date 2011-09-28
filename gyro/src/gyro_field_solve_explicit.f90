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

  call proc_time(CPU_field_in)

  call get_gyro_h
  if (n_field == 3) call get_gyro_h_aperp

  call proc_time(CPU_field2_in)

  ! Velocity-space sums (RHS)
  do ix=1,n_field
     call gyro_velocity_sum(ix)
  enddo

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

  call gyro_field_interpolation

  call proc_time(CPU_field_out)
  CPU_field = CPU_field + (CPU_field_out - CPU_field_in)

end subroutine gyro_field_solve_explicit
