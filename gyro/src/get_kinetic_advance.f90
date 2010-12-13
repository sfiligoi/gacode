!-----------------------------------------------------------
! get_kinetic_advance.f90 [caller: timestep_SSP_322]
!
! PURPOSE:
!  Manage implicit advance of fields and electron
!  distribution.  Explicit ion RHS done with gyro_rhs_total.
!  
! NOTES:
!  Analogous routine for adiabatic electrons is 
!  get_adiabatic_advance. 
!-----------------------------------------------------------

subroutine get_kinetic_advance

  use gyro_globals

  !-------------------------
  implicit none
  !-------------------------

  call proc_time(CPU_field_in)

  call get_delta_he

  call get_gyro_h

  if (n_field == 1) then
     call get_poisson_solution
  else
     call get_maxwell_solution
  endif

  call gyro_field_interpolation

  call get_he

  call proc_time(CPU_field_out)
  CPU_field = CPU_field + (CPU_field_out - CPU_field_in)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[get_kinetic_advance done]'
  endif

end subroutine get_kinetic_advance
