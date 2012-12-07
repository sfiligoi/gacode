!-----------------------------------------------------------
! gyro_timestep_explicit.f90 
!
! PURPOSE:
!  Explicit integrators.
!-----------------------------------------------------------

subroutine gyro_timestep_explicit

  use gyro_globals
  use gyro_pointers

  !--------------------------------------------------------------
  implicit none
  !
  complex, dimension(:,:,:,:), allocatable :: rhse_1
  complex, dimension(:,:,:,:), allocatable :: rhse_2
  complex, dimension(:,:,:,:), allocatable :: rhse_3
  !--------------------------------------------------------------

  allocate(rhse_1(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(rhse_2(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(rhse_3(n_stack,n_x,n_nek_loc_1,n_kinetic))


  if (integrator_method == 1) then

     !---------
     ! RK2
     !---------

     h_0 = h

     ! Stage 1

     call gyro_rhs_total
     rhse_1 = rhs

     h = h_0+dt*rhse_1
     call gyro_field_solve_explicit

     ! Stage 2

     call gyro_rhs_total
     rhse_2 = rhs

     rhs = 0.5*(rhse_1+rhse_2)
     h   = h_0+dt*rhs
     call gyro_field_solve_explicit

     ! Error estimate

     h_err = 0.5*dt*(-rhse_1+rhse_2)

  else 

     !---------
     ! RK4
     !---------

     h_0 = h

     ! Stage 1

     call gyro_rhs_total
     rhse_1 = rhs

     h = h_0+0.5*dt*rhse_1
     call gyro_field_solve_explicit

     ! Stage 2

     call gyro_rhs_total
     rhse_2 = rhs

     h = h_0+0.5*dt*rhse_2
     call gyro_field_solve_explicit

     ! Stage 3

     call gyro_rhs_total
     rhse_3 = rhs

     h = h_0+dt*rhse_3
     call gyro_field_solve_explicit

     ! Stage 4

     call gyro_rhs_total

     h = h_0+dt*(rhse_1+2.0*(rhse_2+rhse_3)+rhs)/6.0
     call gyro_field_solve_explicit

     ! Error estimate

     h_err = dt*(-rhse_1+rhse_2+rhse_3-rhs)/12.0

  endif

  deallocate(rhse_1)
  deallocate(rhse_2)
  deallocate(rhse_3)

end subroutine gyro_timestep_explicit
