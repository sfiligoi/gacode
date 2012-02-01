  !--------------------------------------------------------------
  ! tgyro_global_flux.f90
  !
  ! PURPOSE:
  !  Capture fluxes from global GYRO simulation
  !---------------------------------------------------------------

subroutine tgyro_global_flux

  use tgyro_globals
  use gyro_interface

  implicit none
  
  integer :: i
  integer :: j
  integer :: i_ion
  integer :: n_gyro
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

  do i=2,n_r
     n_gyro = 0
     do j=igmin,igmax
        x = gyro_r_out(j)-r(i)/r_min
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
