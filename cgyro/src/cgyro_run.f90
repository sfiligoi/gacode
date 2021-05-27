!---------------------------------------------------------------
! cgyro_run.f90
!
! PURPOSE:
!  This is the actual system call to run CGYRO via subroutine 
!  interface (i.e, from TGYRO).
!
! NOTES:
! - To run CGYRO as a subroutine, do this:
!   call cgyro_init
!   call cgyro_run
!-------------------------------------------------------------

subroutine cgyro_run(test_flag_in,var_in,n_species_out,flux_tave_out,tave_min_out,tave_max_out)

  use mpi
  use cgyro_globals
  use timer_lib

  implicit none

  ! Input parameters (IN) - REQUIRED
  integer, intent(in)  :: test_flag_in
  integer, intent(in)  :: var_in           
  integer, intent(out) :: n_species_out
  real, intent(out)    :: tave_min_out, tave_max_out
  real, intent(out)    :: flux_tave_out(11,3)

  ! Set corresponding global variables
  test_flag = test_flag_in

  ! Re-set max time
  max_time = var_in
  
  ! Run GYRO
  call cgyro_kernel

  ! Return time-averaged flux data
  n_species_out = n_species
  tave_min_out = tave_min
  tave_max_out = tave_max
  flux_tave_out(:,:) = 0.0
  if(abs(gamma_e) > 1e-10) then
     flux_tave_out(1:n_species,:) = cflux_tave(1:n_species,:)/(1.0*tave_step)
  else
     flux_tave_out(1:n_species,:) = gflux_tave(1:n_species,:)/(1.0*tave_step)
  endif
  
  ! Clean-up
  call cgyro_cleanup
  
end subroutine cgyro_run
