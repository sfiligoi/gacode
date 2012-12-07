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
  complex, dimension(:,:,:,:), allocatable :: rhse_2
  complex, dimension(:,:,:,:), allocatable :: rhse_3
  complex, dimension(:,:,:,:), allocatable :: rhsi_1
  complex, dimension(:,:,:,:), allocatable :: rhsi_2
  complex, dimension(:,:,:,:), allocatable :: rhsi_3
  !--------------------------------------------------------------

  allocate(rhse_2(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(rhse_3(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(rhsi_1(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(rhsi_2(n_stack,n_x,n_nek_loc_1,n_kinetic))
  allocate(rhsi_3(n_stack,n_x,n_nek_loc_1,n_kinetic))


  h_0 = h

  call gyro_kinetic_advance
  rhsi_1 = (h-h_0)/(0.5*dt)

  h = h_0-0.5*dt*rhsi_1

  rhsi_2 = h
  call gyro_kinetic_advance
  rhsi_2 = (h-rhsi_2)/(0.5*dt)

  call gyro_rhs_total
  rhse_2 = rhs

  h = h_0+dt*rhse_2+0.5*dt*rhsi_2

  rhsi_3 = h
  call gyro_kinetic_advance
  rhsi_3 = (h-rhsi_3)/(0.5*dt)

  call gyro_rhs_total
  rhse_3 = rhs

  h = h_0+dt*0.5*(rhse_2+rhse_3+rhsi_2+rhsi_3)
  call gyro_field_solve_explicit

  ! Error estimate

  h_err = 0.5*dt*(-rhse_2+rhse_3-rhsi_2+rhsi_3)

  !----------------------------------------------------------
  ! Use RK4 for nonlinear advance
  !
  if (n_substep > 0 .and. nonlinear_flag == 1) then

     dts = dt/n_substep

     do i_substep=1,n_substep

        h_0 = h
        call gyro_rhs_nonlinear
        rhsi_1 = rhs

        h = h_0+0.5*dts*rhsi_1
        call gyro_field_solve_explicit
        call gyro_rhs_nonlinear
        rhsi_2 = rhs

        h = h_0+0.5*dts*rhsi_2
        call gyro_field_solve_explicit
        call gyro_rhs_nonlinear
        rhsi_3 = rhs

        h = h_0+dts*rhsi_3
        call gyro_field_solve_explicit
        call gyro_rhs_nonlinear

        h = h_0+dts*(rhsi_1+2.0*(rhsi_2+rhsi_3)+rhs)/6.0
        call gyro_field_solve_explicit

        h_err = h_err+dts*(-rhsi_1+rhsi_2+rhsi_3-rhs)/12.0

     enddo

  endif
  !----------------------------------------------------------

  deallocate(rhse_2)
  deallocate(rhse_3)
  deallocate(rhsi_1)
  deallocate(rhsi_2)
  deallocate(rhsi_3)

end subroutine gyro_timestep_implicit
