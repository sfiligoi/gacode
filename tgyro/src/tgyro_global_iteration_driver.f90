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

  integer :: imin,imax
  integer :: j
  integer :: i_ion
  integer :: n_gyro
  real :: length
  real :: dlength
  real :: x

  ! Initialize GYRO
  call gyro_init(paths(1),MPI_COMM_WORLD)

  n_r = tgyro_global_radii+1

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

  tgyro_rmin = gyro_r_out(imin)
  tgyro_rmax = gyro_r_out(imax)
  ! Overwrite TGYRO variables with GYRO values
  length  = tgyro_rmax-tgyro_rmin 
  dlength = length/tgyro_global_radii

  ! Compute internal TGYRO radii
  r(1) = 0.0
  do i=2,n_r
     ! Normalized r (dimensionless)
     r(i) = tgyro_rmin+dlength/2+(i-2)*dlength
  enddo

  print *,tgyro_rmin
  print *,tgyro_rmax

  ! Compute binned fluxes
  pflux_i_tur(:,:) = 0.0
  pflux_e_tur(:) = 0.0
  eflux_i_tur(:,:) = 0.0
  eflux_e_tur(:) = 0.0
  mflux_i_tur(:,:) = 0.0
  mflux_e_tur(:) = 0.0
  expwd_i_tur(:,:) = 0.0
  expwd_e_tur(:) = 0.0
  do i=2,n_r
     n_gyro = 0
     do j=imin,imax
        x = gyro_r_out(j)-r(i)
        ! See if GYRO simulation point is inside TGYRO bin
        if (x > -dlength/2 .and. x < dlength/2) then
           n_gyro = n_gyro+1
           pflux_e_tur(i) = pflux_e_tur(i)+&
                gyro_elec_pflux_out(j)
           eflux_e_tur(i) = eflux_e_tur(i)+&
                gyro_elec_eflux_out(j)  
           mflux_e_tur(i) = mflux_e_tur(i)+&
                gyro_elec_mflux_out(j)  
           expwd_e_tur(i) = expwd_e_tur(i)+&
                gyro_elec_expwd_out(j)  
           do i_ion=1,loc_n_ion
              pflux_i_tur(i_ion,i) = pflux_i_tur(i_ion,i)+&
                   gyro_ion_pflux_out(j,i_ion) 
              eflux_i_tur(i_ion,i) = eflux_i_tur(i_ion,i)+&
                   gyro_ion_eflux_out(j,i_ion)  
              mflux_i_tur(i_ion,i) = mflux_i_tur(i_ion,i)+&
                   gyro_ion_mflux_out(j,i_ion)  
              expwd_i_tur(i_ion,i) = expwd_i_tur(i_ion,i)+&
                   gyro_ion_expwd_out(j,i_ion)  
           enddo ! i_ion
        endif
     enddo ! j
     if (n_gyro == 0) then
        call tgyro_catch_error('ERROR: (TGYRO) TGYRO_GLOBAL_RADII too large.')
     endif
     ! Compute final fluxes by dividing by number of GYRO points
     pflux_e_tur(i) = pflux_e_tur(i)/n_gyro
     eflux_e_tur(i) = eflux_e_tur(i)/n_gyro
     mflux_e_tur(i) = mflux_e_tur(i)/n_gyro
     expwd_e_tur(i) = expwd_e_tur(i)/n_gyro
     do i_ion=1,loc_n_ion
        pflux_i_tur(i_ion,i) = pflux_i_tur(i_ion,i)/n_gyro
        eflux_i_tur(i_ion,i) = eflux_i_tur(i_ion,i)/n_gyro
        mflux_i_tur(i_ion,i) = mflux_i_tur(i_ion,i)/n_gyro
        expwd_i_tur(i_ion,i) = expwd_i_tur(i_ion,i)/n_gyro
     enddo ! i_ion
     print *,i,r(i),eflux_e_tur(i)
  enddo ! i
  print *,sum(eflux_e_tur(2:n_r))/tgyro_global_radii

  call EXPRO_write_original(' ')
  call EXPRO_palloc(MPI_COMM_WORLD,'./',0)
  !--------------------------------------------------------

end subroutine tgyro_global_iteration_driver
