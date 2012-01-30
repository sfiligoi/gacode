!-----------------------------------------------------------
! tgyro_global_iteration_driver.f90
!
! PURPOSE:
!  Main driver for global solver.
!----------------------------------------------------------

subroutine tgyro_global_iteration_driver

  use mpi
  use tgyro_globals
  use tgyro_iteration_variables
  use EXPRO_interface
  use gyro_interface

  implicit none

  integer :: imin,imax,j
  integer :: dn_gyro
  real :: length
  real :: dlength

  ! Initialize GYRO
  call gyro_init(paths(1),MPI_COMM_WORLD)

  n_r = tgyro_global_radii+1

  allocate(f_vec(n_r))
  call tgyro_allocate_globals

  !---------------------------------------
  ! Error monitoring variables
  !
  error_flag = 0
  error_msg = 'INFO: clean exit from TGYRO'
  !
  gyro_exit_status(:)  = 0
  gyro_exit_message(:) = 'N/A'
  !---------------------------------------

  gyro_restart_method = 1
  transport_method    = 2

  call gyro_run(gyrotest_flag, gyro_restart_method, &
       transport_method, gyro_exit_status(1), gyro_exit_message(1))

  !--------------------------------------------------------
  ! Global TGYRO

  EXPRO_ctrl_density_method = loc_quasineutral_flag+1
  EXPRO_ctrl_z = 0.0
  EXPRO_ctrl_z(1:loc_n_ion) = zi_vec(1:loc_n_ion)
  EXPRO_ctrl_numeq_flag = loc_num_equil_flag
  EXPRO_ctrl_signq = tgyro_ipccw_in*tgyro_btccw_in
  EXPRO_ctrl_signb = -tgyro_btccw_in
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_palloc(MPI_COMM_WORLD,'./',1) 
  call EXPRO_pread

  ! GYRO gridpoints corresponding to simulation domain ends
  imin = 1+gyro_explicit_damp_grid_in
  imax = gyro_radial_grid_in-gyro_explicit_damp_grid_in

  ! Overwrite TGYRO variables with GYRO values
  length  = gyro_r_out(imax)-gyro_r_out(imin) 
  dlength = length/tgyro_global_radii
  ! Number of GYRO gridpoints per TGYRO bin
  dn_gyro = (imax-imin)/tgyro_global_radii
  tgyro_rmin = gyro_r_out(imin)
  tgyro_rmax = gyro_r_out(imax)

  ! Compute internal TGYRO radii
  r(1) = 0.0
  do i=2,n_r
     ! Normalized r (dimensionless)
     r(i) = length/2+(i-2)*dlength
  enddo

  ! Compute binned fluxes
  f_vec(:) = 0.0
  do i=2,n_r
     do j=imin+(i-2)*dn_gyro,imin+(i-1)*dn_gyro
        f_vec(i) = f_vec(i)+gyro_elec_eflux_out(i) 
     enddo
     f_vec(i) = f_vec(i)/dn_gyro
     print *,i,f_vec(i)
  enddo

  ! Compute internal TGYRO fluxes
  !  do i=1,gyro_radial_grid_in
  !     print *,gyro_r_out(i),gyro_elec_eflux_out(i)
  !  enddo

  ! Overlay profiles

  ! XXXXX

  call EXPRO_write_original(' ')
  call EXPRO_palloc(MPI_COMM_WORLD,'./',0)
  !--------------------------------------------------------

end subroutine tgyro_global_iteration_driver
