  !--------------------------------------------------------------
  ! tgyro_global_flux.f90
  !
  ! PURPOSE:
  !  Capture fluxes from global GYRO simulation
  !
  ! NOTE:  Exchanges need normalization correction.
  !---------------------------------------------------------------

subroutine tgyro_global_flux

  use tgyro_globals
  use gyro_interface

  implicit none

  integer :: i
  integer :: j
  integer :: i_ion
  integer :: n_gyro
  integer :: inorm
  real :: x

  ! Initialize all fluxes
  pflux_i_neo(:,:) = 0.0
  pflux_e_neo(:) = 0.0
  eflux_i_neo(:,:) = 0.0
  eflux_e_neo(:) = 0.0
  mflux_i_neo(:,:) = 0.0
  mflux_e_neo(:) = 0.0

  pflux_i_tur(:,:) = 0.0
  pflux_e_tur(:) = 0.0
  eflux_i_tur(:,:) = 0.0
  eflux_e_tur(:) = 0.0
  mflux_i_tur(:,:) = 0.0
  mflux_e_tur(:) = 0.0
  expwd_i_tur(:,:) = 0.0
  expwd_e_tur(:) = 0.0

  gyro_restart_method = 1
  call gyro_run(gyrotest_flag, gyro_restart_method, &
       transport_method, gyro_exit_status(1), gyro_exit_message(1))

  !------------------------------------------------------------------
  ! Note that GYRO global fluxes are returned with normalizations  
  ! based on GYRO norm radius, whereas in TGYRO the fluxes are based 
  ! on the GB normalizations at the local radius.  So, we also 
  ! perform this conversion below.
  !
  ! CH: since using an total even number of points
  ! n_r = 1 + tgyro_global_radii want to use

  inorm = n_r/2 + 1

  ! to specify (approximate) center of box.  For instance, with
  ! tgyro_global_radii = 3, corrsponding to [r1,r2,r3], total of
  ! n_r = 4 point (r = [0,r1,r2,r3]) in solver.  Since r2 is norm
  ! point, want inorm = n_r/2 + 1 = 3.
  !------------------------------------------------------------------

  do i=2,n_r
     n_gyro = 0

     do j=igmin,igmax
        x = gyro_r_out(j)-r(i)/r_min
        ! See if GYRO simulation point is inside TGYRO bin
        if (x > -dlength/2 .and. x < dlength/2) then
           n_gyro = n_gyro+1
           pflux_e_tur(i) = pflux_e_tur(i)+&
                gyro_elec_pflux_out(j)*gamma_gb(inorm)/gamma_gb(i)
           eflux_e_tur(i) = eflux_e_tur(i)+&
                gyro_elec_eflux_out(j)*q_gb(inorm)/q_gb(i) 
           mflux_e_tur(i) = mflux_e_tur(i)+&
                gyro_elec_mflux_out(j)*pi_gb(inorm)/pi_gb(i)  
           expwd_e_tur(i) = expwd_e_tur(i)+&
                gyro_elec_expwd_out(j)*0.0
           do i_ion=1,loc_n_ion
              pflux_i_tur(i_ion,i) = pflux_i_tur(i_ion,i)+&
                   gyro_ion_pflux_out(j,i_ion)*gamma_gb(inorm)/gamma_gb(i) 
              eflux_i_tur(i_ion,i) = eflux_i_tur(i_ion,i)+&
                   gyro_ion_eflux_out(j,i_ion)*q_gb(inorm)/q_gb(i) 
              mflux_i_tur(i_ion,i) = mflux_i_tur(i_ion,i)+&
                   gyro_ion_mflux_out(j,i_ion)*pi_gb(inorm)/pi_gb(i) 
              expwd_i_tur(i_ion,i) = expwd_i_tur(i_ion,i)+&
                   gyro_ion_expwd_out(j,i_ion)*0.0  
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
     !print *,i,r(i)/r_min,eflux_e_tur(i)
  enddo ! i
  !print *,sum(eflux_e_tur(2:n_r))/tgyro_global_radii

  !--------------------------------------------------------

  !----------------------------------------------------------
  ! Compute total fluxes given neoclassical and turbulent 
  ! components:
  !
  do i=1,n_r
     pflux_e_tot(i) = pflux_e_neo(i)+pflux_e_tur(i)
     pflux_i_tot(i) = sum(pflux_i_neo(:,i)+pflux_i_tur(:,i))

     eflux_e_tot(i) = eflux_e_neo(i)+eflux_e_tur(i)
     eflux_i_tot(i) = sum(eflux_i_neo(:,i)+eflux_i_tur(:,i))

     mflux_tot(i) = mflux_e_neo(i)+mflux_e_tur(i)+&
          sum(mflux_i_neo(:,i)+mflux_i_tur(:,i))
  enddo
  !----------------------------------------------------------

end subroutine tgyro_global_flux
