!-----------------------------------------------------
! get_field_explicit.f90
!
! PURPOSE:
!  This routine coordinates calling of the appropriate 
!  sequence of routines for solution of the EXPLICIT 
!  Maxwell equations
!-----------------------------------------------------

subroutine get_field_explicit

  use gyro_globals

  implicit none

  call proc_time(CPU_field_in)

  call get_gyro_h

  !---------------------------------------------
  ! Select which field solvers will be used:
  !
  call proc_time(CPU_field2_in)

  if(n_field == 3) then
     call get_maxwell3_explicit
  else
     call get_poisson_explicit
     if (n_field > 1) then
        call get_ampere_explicit
     endif
  endif

  call gyro_field_interpolation
  !---------------------------------------------

  call proc_time(CPU_field_out)
  CPU_field = CPU_field + (CPU_field_out - CPU_field_in)

end subroutine get_field_explicit
