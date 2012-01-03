!-----------------------------------------------------------
! gyro_timestep_implicit.f90 
!
! PURPOSE:
!  High level control of an IMEX RK step.
!
! NOTES:
!  This is the SSP2(3,2,2) IMEX Runge Kutta scheme.
!-----------------------------------------------------------

subroutine gyro_timestep_implicit

  use gyro_globals
  use gyro_pointers

  !--------------------------------------------------------------
  implicit none
  !
  integer :: i_substep
  real :: dts
  !
  complex, dimension(:,:,:,:), allocatable :: RHSE_2
  complex, dimension(:,:,:,:), allocatable :: RHSE_3
  complex, dimension(:,:,:,:), allocatable :: RHSI_1
  complex, dimension(:,:,:,:), allocatable :: RHSI_2
  complex, dimension(:,:,:,:), allocatable :: RHSI_3
  !--------------------------------------------------------------

  allocate(RHSE_2(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(RHSE_3(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(RHSI_1(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(RHSI_2(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(RHSI_3(n_stack,n_x,n_nek_loc_1,n_kinetic))


  h_0 = h

  call get_kinetic_advance
  RHSI_1 = (h-h_0)/(0.5*dt)

  h = h_0-0.5*dt*RHSI_1

  RHSI_2 = h
  call get_kinetic_advance
  RHSI_2 = (h-RHSI_2)/(0.5*dt)

  call gyro_rhs_total
  RHSE_2 = RHS

  h = h_0+dt*RHSE_2+0.5*dt*RHSI_2

  RHSI_3 = h
  call get_kinetic_advance
  RHSI_3 = (h-RHSI_3)/(0.5*dt)

  call gyro_rhs_total
  RHSE_3 = RHS

  h = h_0+dt*0.5*(RHSE_2+RHSE_3+RHSI_2+RHSI_3)
  call gyro_field_solve_explicit

  ! Error estimate

  h_err = 0.5*dt*(-RHSE_2+RHSE_3-RHSI_2+RHSI_3)

  !----------------------------------------------------------
  ! Use RK4 for nonlinear advance
  !
  if (n_substep > 0 .and. nonlinear_flag == 1) then

     dts = dt/n_substep

     do i_substep=1,n_substep

        h_0 = h
        call gyro_rhs_nonlinear
        RHSI_1 = RHS

        h = h_0+0.5*dts*RHSI_1
        call gyro_field_solve_explicit
        call gyro_rhs_nonlinear
        RHSI_2 = RHS

        h = h_0+0.5*dts*RHSI_2
        call gyro_field_solve_explicit
        call gyro_rhs_nonlinear
        RHSI_3 = RHS

        h = h_0+dts*RHSI_3
        call gyro_field_solve_explicit
        call gyro_rhs_nonlinear

        h = h_0+dts*(RHSI_1+2.0*(RHSI_2+RHSI_3)+RHS)/6.0
        call gyro_field_solve_explicit

        h_err = h_err+dts*(-RHSI_1+RHSI_2+RHSI_3-RHS)/12.0

     enddo

  endif
  !----------------------------------------------------------

  deallocate(RHSE_2)
  deallocate(RHSE_3)
  deallocate(RHSI_1)
  deallocate(RHSI_2)
  deallocate(RHSI_3)

end subroutine gyro_timestep_implicit
