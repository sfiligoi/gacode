!---------------------------------------------------------------
! gyro_run.f90
!
! PURPOSE:
!  This is the actual system call to run GYRO via subroutine 
!  interface (i.e, from TGYRO).
!
! NOTES:
! - To run GYRO as a subroutine, do this:
!   call gyro_init
!   call gyro_run
!  
! - To run inside TGYRO, we do this:
!   call tgyro_gyro_map (calls gyro_init)
!   call gyro_run(...)
!-------------------------------------------------------------

subroutine gyro_run(&
     test_flag_in,&
     restart_method_in,&
     transport_method_in,&
     status_out,&
     message_out)

  use mpi
  use gyro_globals
  use gyro_interface
  use gyro_fieldeigen_private

  implicit none

  ! Input parameters (IN) - REQUIRED
  integer, intent(in) :: test_flag_in
  integer, intent(in) :: restart_method_in
  integer, intent(in) :: transport_method_in

  ! Input parameters (OUT) - REQUIRED
  integer,          intent(inout) :: status_out
  character(len=*), intent(inout) :: message_out

  ! Local variables
  integer :: i_ion
  integer :: err

  ! Set corresponding global variables
  gyrotest_flag    = test_flag_in
  restart_method   = restart_method_in
  transport_method = transport_method_in

  ! NOTE:
  !   (1) ion temperature information is currently automatically equated below
  !       only if 'undefined' (referring to DLNTDR_2=3.0, DLNTDR_3=3.0,
  !       TI_OVER_TE_2=1.0, and TI_OVER_TE_3=1.0, respectively)
  !   (2) ion species temperature information can be explicitly set by setting
  !       the values for the appropriate variables thus overriding the default
  !       behavior

  ! Set Tz = Ti for all ion species
  if (gyro_ni_over_ne_2_in == 0.0) then
     gyro_dlntdr_2_in = gyro_dlntdr_in
     gyro_ti_over_te_2_in = gyro_ti_over_te_in
  endif
  if (gyro_ni_over_ne_3_in == 0.0) then
     gyro_dlntdr_3_in = gyro_dlntdr_in
     gyro_ti_over_te_3_in = gyro_ti_over_te_in
  endif

  ! Map INTERFACE parameters -> GLOBAL variables and allocate output arrays
  call map_interface2global

  ! Run GYRO
  call gyro_do
  call send_line('STATUS: '//gyro_exit_message)

  ! Get output information

  status_out  = gyro_exit_status
  message_out = gyro_exit_message

  ! Running in test mode
  if (gyrotest_flag == 1) then
     return
  endif

  ! GYRO failed
  if (gyro_exit_status == 1) then
     return
  endif

  ! Linear stability (FIELDEIGEN)
  if (linsolve_method == 3) then

     gyro_fieldeigen_omega_out = omega_eigen
     gyro_fieldeigen_error_out = error_eigen

     call gyro_cleanup
     return

  endif

  !-------------------------------------------------------------------------------------- 
  ! Pack transport fluxes into interface output arrays
  !
  if (i_proc == 0 .and. transport_method > 1) then

     if (electron_method == 2) then
        do i=1,n_x
           ! Gamma_e/Gamma_GB
           gyro_elec_pflux_out(i) = sum(gbflux_vec(indx_e,1:n_field,1,i))
           ! Q_e/Q_GB
           gyro_elec_eflux_out(i) = sum(gbflux_vec(indx_e,1:n_field,2,i))
           ! Pi_e/Pi_GB
           gyro_elec_mflux_out(i) = sum(gbflux_vec(indx_e,1:n_field,3,i))
           ! S_e/S_GB
           gyro_elec_expwd_out(i) = sum(gbflux_vec(indx_e,1:n_field,4,i))
        enddo
     endif

     do i_ion=1,n_ion
        do i=1,n_x
           ! Gamma_i/Gamma_GB
           gyro_ion_pflux_out(i,i_ion) = sum(gbflux_vec(i_ion,1:n_field,1,i))
           ! Q_i/Q_GB
           gyro_ion_eflux_out(i,i_ion) = sum(gbflux_vec(i_ion,1:n_field,2,i))
           ! Pi_i/Pi_GB
           gyro_ion_mflux_out(i,i_ion) = sum(gbflux_vec(i_ion,1:n_field,3,i))
           ! S_i/S_GB
           gyro_ion_expwd_out(i,i_ion) = sum(gbflux_vec(i_ion,1:n_field,4,i))
        enddo
     enddo

     gyro_r_out(:) = r(:)

  endif
  !-------------------------------------------------------------------------------------- 

  !-------------------------------------------------------------------------------------- 
  ! Broadcast results to all cores:
  !
  ! Electrons
  call MPI_BCAST(gyro_elec_pflux_out,n_x,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)
  call MPI_BCAST(gyro_elec_eflux_out,n_x,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)
  call MPI_BCAST(gyro_elec_mflux_out,n_x,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)
  call MPI_BCAST(gyro_elec_expwd_out,n_x,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)

  ! Ions
  call MPI_BCAST(gyro_ion_pflux_out,n_x*n_ion,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)
  call MPI_BCAST(gyro_ion_eflux_out,n_x*n_ion,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)
  call MPI_BCAST(gyro_ion_mflux_out,n_x*n_ion,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)
  call MPI_BCAST(gyro_ion_expwd_out,n_x*n_ion,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)

  call MPI_BCAST(gyro_r_out,n_x,MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,err)
  !-------------------------------------------------------------------------------------- 

  call gyro_cleanup

  if (debug_flag == 1 .and. i_proc == 0) then
     print *, '[gyro_run done]'
  endif

end subroutine gyro_run
