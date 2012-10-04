!---------------------------------------------------------------------
! tgyro_global_iteration_driver.f90-
!
! PURPOSE:
!  Main driver for global solver.
!
! NOTES:
!  Management of input.profiles is somewhat complicated.  We will end 
!  up with the following structure:
!
!  tgyrodir/input.profiles
!  tgyrodir/input.profiles.gen
!  tgyrodir/input.profiles.{n}
!  tgyrodir/input.profiles.{n}.gen
!
!  tgyrodir/GYRO*/input.profiles.gen
!  
!  In general, the GYRO*/input.profiles.gen will correspond to the 
!  latest tgyrodir/input.profiles.{n}.gen.
!
!  GLOBAL-TO-LOCAL GRID MAPPING:
!
!  n_r = 6
!  tgyro_global_radii = 5
!
!                            length
!                  <-------------------------> 
!  |
!  |            left GYRO                right GYRO
!  |             buffer                    buffer
!  |               +                         +
!  |               +            dlength      +
!  |               +            <---->       +
!  |               +  |    |    |    |    |  +
!  ------------------------------------------------------
! i=1                i=2  i=3  i=4  i=5  i=6
! r=0                
!
!---------------------------------------------------------------------

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

  p_max = n_evolve*(n_r-1)

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

  ! Map function from radius/field to p-index.
  p = 0
  do ip=1,n_evolve
     do i=2,n_r
        p = p+1
        pmap(i,ip) = p
     enddo
  enddo

  ! Set initial values
  res(:)   = 1.0
  relax(:) = 0.0

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
  !  i_tran initialized in tgyro_global_init_profiles
  !  If not restarting, reduce by 1 to get iteration 0 in loop
  !
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

     i_tran = i_tran+1

     ! Integrate profiles based on gradients
     call tgyro_profile_functions

     ! Read profile data, copy current profiles into interface, 
     ! rewrite profiles
     call EXPRO_palloc(MPI_COMM_WORLD,'./',1) 
     call EXPRO_pread

     !------------------------------------------------------------------------------------------
     ! TGYRO-TO-GYRO PROFILE MAPPING:
     !
     ! Map Te,ze from TGYRO variable to EXPRO interface variable
     call tgyro_global_interpolation(r/1e2,dlntedr*1e2,te/1e3,&
          n_r,n_exp,EXPRO_rmin,EXPRO_te)

     ! Map Ti,zi from TGYRO variable to EXPRO interface variable
     call tgyro_global_interpolation(r/1e2,dlntidr(1,:)*1e2,ti(1,:)/1e3,& 
          n_r,n_exp,EXPRO_rmin,EXPRO_ti(1,:))

     EXPRO_ti(2,:) = EXPRO_ti(1,:)

     ! Map ne,zne from TGYRO variable to EXPRO interface variable
     call tgyro_global_interpolation(r/1e2,dlnnedr*1e2,ne/1e13,n_r,&
          n_exp,EXPRO_rmin,EXPRO_ne)

     ! Map ni,zni from TGYRO variable to EXPRO interface variable, enforcing quasineutrality as appropriate
     if (loc_quasineutral_flag == 1) then
        call tgyro_quasineutral(ni,ne,dlnnidr,dlnnedr,zi_vec,loc_n_ion,n_r)
     endif
     call tgyro_global_interpolation(r/1e2,dlnnidr(1,:)*1e2,ni(1,:)/1e13,n_r,&
          n_exp,EXPRO_rmin,EXPRO_ni(1,:))
     !------------------------------------------------------------------------------------------

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
     ! 
     !  (delta z)/z = - alpha (delta Q)/Q
     !  z = (1 - alpha (dQ/Q))*z
     ! 
     ! alpha = 1/(Waltz stiffness), use alpha=0.1 <=> S = 10
     ! use alpha = loc_relax as TGYRO input

     p = 0

     if (loc_ti_feedback_flag == 1) then
        do i=2,n_r
           p = p+1
           res(p) = (eflux_i_tot(i)-eflux_i_target(i))/eflux_i_target(i) 
           dlntidr(1,i) = dlntidr(1,i)-loc_relax*res(p)*dlntidr(1,i)
        enddo
     endif

     if (loc_te_feedback_flag == 1) then
        do i=2,n_r
           p = p+1
           res(p) = (eflux_e_tot(i)-eflux_e_target(i))/eflux_e_target(i)
           dlntedr(i) = dlntedr(i)-loc_relax*res(p)*dlntedr(i)
        enddo
     endif

     if (loc_ne_feedback_flag == 1) then
        do i=2,n_r
           p = p+1
           res(p) = (pflux_e_tot(i)-pflux_e_target(i))/max(abs(pflux_e_tot(i)),1.0) 
           dlnnedr(i) = dlnnedr(i)-loc_relax*res(p)*dlnnedr(i)
        enddo
        if (loc_quasineutral_flag == 1) then
           call tgyro_quasineutral(ni,ne,dlnnidr,dlnnedr,zi_vec,loc_n_ion,n_r)
        endif
     endif
     !------------------------------------------------------------

     res = abs(res)

  enddo

end subroutine tgyro_global_iteration_driver
