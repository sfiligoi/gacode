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
  complex, dimension(:,:,:,:), allocatable :: RHSE_1
  complex, dimension(:,:,:,:), allocatable :: RHSE_2
  complex, dimension(:,:,:,:), allocatable :: RHSE_3
  !--------------------------------------------------------------

  allocate(RHSE_1(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(RHSE_2(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(RHSE_3(n_stack,n_x,n_nek_loc_1,n_kinetic))


  if (integrator_method == 1) then

     !---------
     ! RK2
     !---------

     h_0 = h

     ! Stage 1

     call gyro_rhs_total
     RHSE_1 = RHS

     h = h_0+dt*RHSE_1
     call get_field_explicit

     ! Stage 2

     call gyro_rhs_total
     RHSE_2 = RHS

     RHS = 0.5*(RHSE_1+RHSE_2)
     h = h_0+dt*RHS
     call get_field_explicit

     ! Error estimate

     h_err = 0.5*dt*(-RHSE_1+RHSE_2)

  else 

     !---------
     ! RK4
     !---------

     h_0 = h

     ! Stage 1

     call gyro_rhs_total
     RHSE_1 = RHS

     h = h_0+0.5*dt*RHSE_1
     call get_field_explicit

     ! Stage 2

     call gyro_rhs_total
     RHSE_2 = RHS

     h = h_0+0.5*dt*RHSE_2
     call get_field_explicit

     ! Stage 3

     call gyro_rhs_total
     RHSE_3 = RHS

     h = h_0+dt*RHSE_3
     call get_field_explicit

     ! Stage 4

     call gyro_rhs_total

     h = h_0+dt*(RHSE_1+2.0*(RHSE_2+RHSE_3)+RHS)/6.0
     call get_field_explicit

     ! Error estimate

     h_err = dt*(-RHSE_1+RHSE_2+RHSE_3-RHS)/12.0

  endif

  deallocate(RHSE_1)
  deallocate(RHSE_2)
  deallocate(RHSE_3)

end subroutine gyro_timestep_explicit
