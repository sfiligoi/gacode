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
  character (len=16) :: ittag
! CH: replace with loc_relax
!  real, parameter :: cgrad = 1e-20
!  real, parameter :: cgrad = 0.1

  ! Copy (TGYRO copy of input.profiles) -> (GYRO copy of input.profiles)
  if (i_proc_global == 0) then
     call system('$GACODE_ROOT/tgyro/bin/tgyro_global_helper input.profiles '//trim(paths(1)))
  endif

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

  !--------------------------------------------------------
  ! INITIALIZE GLOBAL TGYRO PROFILES
  !
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

  EXPRO_ctrl_density_method = loc_quasineutral_flag+1
  EXPRO_ctrl_z = 0.0
  EXPRO_ctrl_z(1:loc_n_ion) = zi_vec(1:loc_n_ion)
  EXPRO_ctrl_numeq_flag = loc_num_equil_flag
  EXPRO_ctrl_signq = tgyro_ipccw_in*tgyro_btccw_in
  EXPRO_ctrl_signb = -tgyro_btccw_in
  EXPRO_ctrl_rotation_method = 1

  call EXPRO_palloc(MPI_COMM_WORLD,'./',1) 
  call EXPRO_pread

  n_exp = EXPRO_n_exp

  call tgyro_global_init_profiles

  call EXPRO_write_original('REWROTE')
  call EXPRO_palloc(MPI_COMM_WORLD,'./',0)

  ! Output initialization
! CH i_tran initialized in tgyro_global_init_profiles
!  i_tran = 0
  if (i_tran == 0) then
	i_tran = -1
  endif
  call tgyro_write_input
  call tgyro_write_data(0)
  !--------------------------------------------------------

  !------------------------------------------------------------
  ! TGYRO-GYRO ITERATION CYCLE
  !
  do i_tran_loop=1,tgyro_relax_iterations

     !CH: account for iteration 0 off initial profiles
!     if ((i_tran_loop > 1) then
	i_tran = i_tran+1
!     endif

     ! Integrate profiles based on gradients
     call tgyro_profile_functions

     ! Read profile data, copy current profiles into interface, 
     ! rewrite profiles
     call EXPRO_palloc(MPI_COMM_WORLD,'./',1) 
     call EXPRO_pread

     ! Map Te,ze from TGYRO variable to EXPRO interface variable
     call tgyro_global_interpolation(r/1e2,dlntedr*1e2,te/1e3,n_r,n_exp,EXPRO_rmin,EXPRO_te)

     ! Map Ti,zi from TGYRO variable to EXPRO interface variable
     call tgyro_global_interpolation(r/1e2,dlntidr(1,:)*1e2,ti(1,:)/1e3,n_r,n_exp,EXPRO_rmin,EXPRO_ti(1,:))

     EXPRO_ti(2,:) = EXPRO_ti(1,:)

     ittag = '.'//achar(i_tran_loop-1+iachar("1"))  

     call EXPRO_write_original('REWRITE'//ittag)
     call EXPRO_palloc(MPI_COMM_WORLD,'./',0) 
     if (i_proc_global == 0) then
        call system('cp input.profiles.new input.profiles'//ittag)
        call system('$GACODE_ROOT/tgyro/bin/tgyro_global_helper input.profiles'//ittag//' '//trim(paths(1)))
     endif

     ! Get global GYRO flux, compute targets, write data
     call tgyro_global_flux
     call tgyro_source
     call tgyro_write_data(1)

     !------------------------------------------------------------
     ! MODIFY GRADIENTS
     !
     ! Modify gradient profile based on some "diagonal rule"
     !  dlntedr(:) = 0.0*(eflux_e_tot(:)-eflux_e_target(:))+dlntedr(:)
     ! CH test: try (delta z)/z = - alpha (delta Q)/Q
     !  --> z = (1 - alpha (dQ/Q))*z
     ! alpha = 1/(Waltz stiffness), use alpha=0.1 <=> S = 10
     ! use alpha = loc_relax as TGYRO input

     if (loc_te_feedback_flag == 1) then
     	dlntedr(2:n_r) = -loc_relax*( (eflux_e_tot(2:n_r) &
             - eflux_e_target(2:n_r))/eflux_e_target(2:n_r) )*dlntedr(2:n_r) + dlntedr(2:n_r)
     endif	

     if (loc_te_feedback_flag == 1) then
        dlntidr(1,2:n_r) = -loc_relax*( (eflux_i_tot(2:n_r) &
             - eflux_i_target(2:n_r))/eflux_i_target(2:n_r) )*dlntidr(1,2:n_r) + dlntidr(1,2:n_r)
     endif

     ! Not needed
     !     dlntedr(1)   = 0.0
     !     dlntidr(:,1) = 0.0
     !------------------------------------------------------------

  enddo

end subroutine tgyro_global_iteration_driver
