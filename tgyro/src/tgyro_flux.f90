!------------------------------------------------------------
! tgyro_flux.f90
!
! PURPOSE:
!  Manage calls to obtain neoclassical and turbulent fluxes.
!------------------------------------------------------------

subroutine tgyro_flux

  use mpi
  use tgyro_globals
  use gyro_interface
  use neo_interface
  use tglf_interface
  use qfm_interface

  implicit none

  integer :: i_ion
  integer :: i1,i2
  integer :: n_12
  real :: dummy1
  real :: dummy2
  real :: rltcrit
  real :: rltcritz
  real :: chii0
  real :: chie0
  real :: Gamma_neo_GB
  real :: Q_neo_GB
  real :: Pi_neo_GB
  real, dimension(8) :: x_out

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,'(t2,a)') 'INFO: (TGYRO) Entered tgyro_flux'
     close(1)
  endif

  !-----------------------------
  ! Counter
  flux_counter = flux_counter+1
  !-----------------------------

  !-----------------------------------------------------------
  ! Neoclassical fluxes
  !-----------------------------------------------------------

  pflux_i_neo(:,:) = 0.0
  pflux_e_neo(:)   = 0.0
  eflux_i_neo(:,:) = 0.0
  eflux_e_neo(:)   = 0.0
  mflux_i_neo(:,:) = 0.0
  mflux_e_neo(:)   = 0.0

  ! Map TGYRO parameters to NEO
  call tgyro_neo_map  

  ! NEO normalization to GB output conversions
  ! Gamma_norm = n0_norm vt_norm
  ! Q_norm     = n0_norm vt_norm T_norm
  ! Pi_norm    = n0_norm T0_norm a_norm
  ! where
  ! vt_norm = sqrt(T_norm / mass_norm)
  !
  ! Gamma_GB = n_0e*c_s*(rho_s/a)^2
  ! Q_GB     = n_0e*c_s*(rho_s/a)^2*T_0e
  ! Pi_GB    = n_0e*T_0e*a*(rho_s/a)^2
  ! where
  ! c_s    = sqrt(T_0e/m_i)
  ! rho_se = c_s/Omega_ci
  Gamma_neo_GB = neo_dens_2_in * neo_temp_2_in**1.5 * neo_rho_star_in**2
  Q_neo_GB     = neo_dens_2_in * neo_temp_2_in**2.5 * neo_rho_star_in**2
  Pi_neo_GB    = neo_dens_2_in * neo_temp_2_in**2   * neo_rho_star_in**2

  select case (loc_neo_method) 

  case (1)

     ! analytic theory
     call neo_run

     pflux_i_neo(1,i_r) = neo_pflux_thHH_out/Gamma_neo_GB
     if (loc_chang_hinton == 1) then
        eflux_i_neo(1,i_r) = neo_eflux_thCHi_out/Q_neo_GB
     else
        eflux_i_neo(1,i_r) = neo_eflux_thHHi_out/Q_neo_GB
     endif

     pflux_e_neo(i_r) = neo_pflux_thHH_out /Gamma_neo_GB 
     eflux_e_neo(i_r) = neo_eflux_thHHe_out /Q_neo_GB

  case (2)

     ! dke solve
     if (gyrotest_flag == 0) call neo_run

     pflux_e_neo(i_r)   = (neo_pflux_dke_out(2) + neo_gv_flag*neo_pflux_gv_out(2)) &
          / Gamma_neo_GB 
     eflux_e_neo(i_r)   = (neo_efluxncv_dke_out(2) + neo_gv_flag*neo_efluxncv_gv_out(2)) &
          / Q_neo_GB
     mflux_e_neo(i_r)   = (neo_mflux_dke_out(2) + neo_gv_flag*neo_mflux_gv_out(2)) &
          / Pi_neo_GB

     pflux_i_neo(1,i_r) = (neo_pflux_dke_out(1) + neo_gv_flag*neo_pflux_gv_out(1)) &
          / Gamma_neo_GB 
     eflux_i_neo(1,i_r) = (neo_efluxncv_dke_out(1) + neo_gv_flag*neo_efluxncv_gv_out(1)) &
          / Q_neo_GB
     mflux_i_neo(1,i_r) = (neo_mflux_dke_out(1) + neo_gv_flag*neo_mflux_gv_out(1)) &
          / Pi_neo_GB

     do i_ion=2,loc_n_ion
        pflux_i_neo(i_ion,i_r) = (neo_pflux_dke_out(i_ion+1) + neo_gv_flag*neo_pflux_gv_out(i_ion+1)) &
             / Gamma_neo_GB 
        eflux_i_neo(i_ion,i_r) = (neo_efluxncv_dke_out(i_ion+1) + neo_gv_flag*neo_efluxncv_gv_out(i_ion+1)) &
             / Q_neo_GB
        mflux_i_neo(i_ion,i_r) = (neo_mflux_dke_out(i_ion+1) + neo_gv_flag*neo_mflux_gv_out(i_ion+1)) &
             / Pi_neo_GB
     enddo

  end select

  !-----------------------------------------------------------
  ! Turbulent fluxes
  !-----------------------------------------------------------

  pflux_i_tur(:,:) = 0.0
  pflux_e_tur(:)   = 0.0
  eflux_i_tur(:,:) = 0.0
  eflux_e_tur(:)   = 0.0
  mflux_i_tur(:,:) = 0.0
  mflux_e_tur(:)   = 0.0
  expwd_i_tur(:,:) = 0.0
  expwd_e_tur(:)   = 0.0

  select case (flux_method)

  case (1) 

     call ifs_pppl(r_maj(i_r)*dlntidr(1,i_r),&
          r_maj(i_r)*dlnnidr(1,i_r),&
          r_maj(i_r)*dlnnedr(i_r),&
          abs(q(i_r)),&
          kappa(i_r),&
          s(i_r),&
          1.0,&
          0.0,&
          ti(1,i_r)/te(i_r),&
          r(i_r)/r_maj(i_r),&
          2.5e-7*ne(i_r)/(te(i_r)**1.5*ti(1,i_r)**0.5)*(r_maj(i_r)/100), &
          r_maj(i_r),&
          rho_i(i_r),&
          v_i(i_r), &
          rltcrit,&
          rltcritz,&
          dummy1,&
          dummy2, &
          chii0,&
          chie0)

     ! Normalize to chi_GB = rho_s^2 c_s/a
     x_out(2) = chie0/chi_gb(i_r)
     x_out(4) = chii0/chi_gb(i_r)

     ! Convert to flux
     eflux_e_tur(i_r) = x_out(2)*r_min*dlntedr(i_r)

     eflux_i_tur(1,i_r) = x_out(4)*r_min*dlntidr(1,i_r)*&
          ni(1,i_r)/ne(i_r)*ti(1,i_r)/te(i_r)

  case (2)

     ! Map TGYRO parameters to TGLF
     call tgyro_tglf_map

     if (gyrotest_flag == 0) call tglf_run

     pflux_e_tur(i_r) = tglf_elec_pflux_out
     eflux_e_tur(i_r) = tglf_elec_eflux_out
     mflux_e_tur(i_r) = tglf_elec_mflux_out
     expwd_e_tur(i_r) = tglf_elec_expwd_out

     pflux_i_tur(1:loc_n_ion,i_r) = tglf_ion_pflux_out(1:loc_n_ion)
     eflux_i_tur(1:loc_n_ion,i_r) = tglf_ion_eflux_out(1:loc_n_ion)
     mflux_i_tur(1:loc_n_ion,i_r) = tglf_ion_mflux_out(1:loc_n_ion)
     expwd_i_tur(1:loc_n_ion,i_r) = tglf_ion_expwd_out(1:loc_n_ion)

     if (tglf_q_low_flag == 1) then
        eflux_e_tur(i_r) = tglf_elec_eflux_low_out
        eflux_i_tur(1:loc_n_ion,i_r) = tglf_ion_eflux_low_out(1:loc_n_ion)
     endif

  case (3)

     ! Map TGYRO parameters to GYRO
     call tgyro_gyro_map

     call gyro_run(gyrotest_flag, gyro_restart_method, &
          transport_method, gyro_exit_status(i_r), gyro_exit_message(i_r))

     if (tgyro_mode == 1) then

        i1 = 1+gyro_explicit_damp_grid_in
        i2 = gyro_radial_grid_in-gyro_explicit_damp_grid_in
        n_12 = gyro_radial_grid_in-2*gyro_explicit_damp_grid_in

        ! Map GYRO (local simulation) output to TGYRO

        pflux_e_tur(i_r) = sum(gyro_elec_pflux_out(i1:i2))/n_12
        eflux_e_tur(i_r) = sum(gyro_elec_eflux_out(i1:i2))/n_12
        mflux_e_tur(i_r) = sum(gyro_elec_mflux_out(i1:i2))/n_12
        expwd_e_tur(i_r) = sum(gyro_elec_expwd_out(i1:i2))/n_12

        do i_ion=1,loc_n_ion
           pflux_i_tur(i_ion,i_r) = sum(gyro_ion_pflux_out(i1:i2,i_ion))/n_12
           eflux_i_tur(i_ion,i_r) = sum(gyro_ion_eflux_out(i1:i2,i_ion))/n_12
           mflux_i_tur(i_ion,i_r) = sum(gyro_ion_mflux_out(i1:i2,i_ion))/n_12
           expwd_i_tur(i_ion,i_r) = sum(gyro_ion_expwd_out(i1:i2,i_ion))/n_12
        enddo

     endif

  case default

     call tgyro_catch_error('ERROR: no matching flux method.')

  end select

  !----------------------------------------------------------
  ! Compute total fluxes given neoclassical and turbulent 
  ! components:
  !
  pflux_e_tot(i_r) = pflux_e_neo(i_r)+pflux_e_tur(i_r)
  pflux_i_tot(i_r) = sum(pflux_i_neo(:,i_r)+pflux_i_tur(:,i_r))

  eflux_e_tot(i_r) = eflux_e_neo(i_r)+eflux_e_tur(i_r)
  eflux_i_tot(i_r) = sum(eflux_i_neo(:,i_r)+eflux_i_tur(:,i_r))

  mflux_tot(i_r) = mflux_e_neo(i_r)+mflux_e_tur(i_r)+&
       sum(mflux_i_neo(:,i_r)+mflux_i_tur(:,i_r))
  !----------------------------------------------------------

  !----------------------------------------------------------
  ! Some nonsensical "fudges" to account for sawteeth.  
  ! In general, these should NOT be used.
  !
  select case (loc_sawtooth_model)

  case (2)

     if (q(i_r) <= 1.0) then
        eflux_e_tot(i_r) = eflux_i_neo(1,i_r)+eflux_e_tur(i_r)
     endif

  case (3)

     if (q(i_r) <= 1.0) then
        eflux_e_tot(i_r) = eflux_i_neo(1,i_r)+eflux_e_tur(i_r)
        eflux_i_tot(i_r) = 10*eflux_i_neo(1,i_r)+eflux_i_tur(1,i_r)
     endif

  end select
  !----------------------------------------------------------

end subroutine tgyro_flux
