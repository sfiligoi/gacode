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

  real :: time_max_save
  integer :: n_exp

  ! Initialize GYRO
  call gyro_init(paths(1),MPI_COMM_WORLD)

  n_r = tgyro_global_radii+1
  if (2*int(tgyro_global_radii/2) == tgyro_global_radii) then
     call tgyro_catch_error('ERROR: Must have odd number of GYRO radii')
  endif

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

  ! gyro_restart_method = 0 (no restart)
  !                     = 1 (standard restart)

  gyro_restart_method = 0
  transport_method    = 2
  time_max_save = gyro_time_max_in
  gyro_time_max_in    = 0.0
  call gyro_run(gyrotest_flag, gyro_restart_method, &
       transport_method, gyro_exit_status(1), gyro_exit_message(1))
  gyro_time_max_in    = time_max_save

  ! GYRO gridpoints corresponding to simulation domain ends
  igmin = 1+gyro_explicit_damp_grid_in
  igmax = gyro_radial_grid_in-gyro_explicit_damp_grid_in

  tgyro_rmin = gyro_r_out(igmin)
  tgyro_rmax = gyro_r_out(igmax)
  ! Overwrite TGYRO variables with GYRO values
  length  = tgyro_rmax-tgyro_rmin 
  dlength = length/tgyro_global_radii

  ! Compute internal TGYRO radii
  r(1) = 0.0
  do i=2,n_r
     ! Normalized r (dimensionless)
     r(i) = tgyro_rmin+dlength/2+(i-2)*dlength
  enddo

  !--------------------------------------------------------
  EXPRO_ctrl_density_method = loc_quasineutral_flag+1
  EXPRO_ctrl_z = 0.0
  EXPRO_ctrl_z(1:loc_n_ion) = zi_vec(1:loc_n_ion)
  EXPRO_ctrl_numeq_flag = loc_num_equil_flag
  EXPRO_ctrl_signq = tgyro_ipccw_in*tgyro_btccw_in
  EXPRO_ctrl_signb = -tgyro_btccw_in
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_palloc(MPI_COMM_WORLD,paths(1),1) 
  call EXPRO_pread

  n_exp = EXPRO_n_exp

  call tgyro_global_init_profiles

  call EXPRO_write_original('REWROTE')
  call EXPRO_palloc(MPI_COMM_WORLD,paths(1),0)
  !--------------------------------------------------------

  ! Output initialization
  call tgyro_write_input
  call tgyro_write_data(0)

  ! Integrate profiles based on gradients
  call tgyro_profile_functions

  !------------------------------------------------------------
  ! Read profile data, copy current profiles into interface, 
  ! rewrite profiles
  !
  call EXPRO_palloc(MPI_COMM_WORLD,paths(1),1) 
  call EXPRO_pread

  call cub_spline(r/r_min,te/1e3,n_r,100*EXPRO_rmin(:)/r_min,EXPRO_te(:),n_exp)
  print *,r/r_min,te/1e3
  print *,100*EXPRO_rmin(:)/r_min,EXPRO_te(:)

  call EXPRO_write_original('REWROTE_1')
  call EXPRO_palloc(MPI_COMM_WORLD,paths(1),0)
  !------------------------------------------------------------

  call tgyro_catch_error('ks')

  ! 1: Get global GYRO flux, compute targets, write data
  call tgyro_global_flux
  call tgyro_source
  call tgyro_write_data(1)

  ! Modify gradient profile based on some "diagonal rule"
  !dlntedr(:) = 0.1*(eflux_e_tot(:)-eflux_e_target(:))+dlntedr(:)

  ! Integrate profiles based on gradients
  call tgyro_profile_functions

  !--------------------------------------------------------
  ! Read profile data, copy current profiles into interface, 
  ! rewrite profiles
  !
  call EXPRO_palloc(MPI_COMM_WORLD,paths(1),1) 
  call EXPRO_pread

  call cub_spline(r,te/1e3,n_r,EXPRO_rmin(:)/r_min,EXPRO_te(:),n_exp)

  call EXPRO_write_original('REWROTE_1')
  call EXPRO_palloc(MPI_COMM_WORLD,paths(1),0)
  !--------------------------------------------------------

  ! 2: Get global GYRO flux, compute targets, write data
  call tgyro_global_flux
  call tgyro_source
  call tgyro_write_data(1)

end subroutine tgyro_global_iteration_driver
