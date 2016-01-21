!---------------------------------------------------------
! tglf_run_mpi.f90
!
! PURPOSE:
!  Manage call to TGLF simulation for both standalone and
!  TGYRO usage.
!---------------------------------------------------------

subroutine tglf_run_nn()

  use tglf_pkg
  use tglf_interface
  use tglf_mpi

  implicit none

  real :: OUT_ENERGY_FLUX_1, OUT_PARTICLE_FLUX_1, OUT_STRESS_TOR_1
  real :: OUT_ENERGY_FLUX_3, OUT_PARTICLE_FLUX_3, OUT_STRESS_TOR_3

  integer :: i_ion,n
  complex :: xi=(0.0,1.0)

  call put_signs(tglf_sign_Bt_in,tglf_sign_It_in)

  call put_species(tglf_ns_in, &
       tglf_zs_in, &
       tglf_mass_in)

  call put_kys(tglf_ky_in)

  call put_gaussian_width(tglf_width_in, &
       tglf_width_min_in, &
       tglf_nwidth_in, &
       tglf_find_width_in)

  call put_eikonal(tglf_new_eikonal_in)

  call put_gradients(tglf_rlns_in, &
       tglf_rlts_in, &
       tglf_vpar_shear_in, &
       tglf_vexb_shear_in)

  call put_averages(tglf_taus_in, &
       tglf_as_in, &
       tglf_vpar_in, &
       tglf_vexb_in, &
       tglf_betae_in, &
       tglf_xnue_in, &
       tglf_zeff_in, &
       tglf_debye_in)

  call put_switches(tglf_iflux_in, &
       tglf_use_bper_in, &
       tglf_use_bpar_in, &
       tglf_use_mhd_rule_in, &
       tglf_use_bisection_in, &
       tglf_ibranch_in, &
       tglf_nmodes_in, &
       tglf_nbasis_max_in, &
       tglf_nbasis_min_in, &
       tglf_nxgrid_in, &
       tglf_nky_in)

  call put_rare_switches(tglf_theta_trapped_in, &
       tglf_park_in, &
       tglf_ghat_in, &
       tglf_gchat_in, &
       tglf_wd_zero_in, &
       tglf_linsker_factor_in, &
       tglf_gradB_factor_in, &
       tglf_filter_in, &
       tglf_damp_psi_in, &
       tglf_damp_sig_in)

  call put_model_parameters(tglf_adiabatic_elec_in, &
       tglf_alpha_e_in, &
       tglf_alpha_p_in, &
       tglf_alpha_mach_in, &
       tglf_alpha_quench_in, &
       tglf_alpha_zf_in, &
       tglf_xnu_factor_in, &
       tglf_debye_factor_in, &
       tglf_etg_factor_in, &
       tglf_sat_rule_in, &
       tglf_kygrid_model_in, &
       tglf_xnu_model_in, &
       tglf_vpar_model_in, &
       tglf_vpar_shear_model_in)

  if (tglf_geometry_flag_in == 1 ) then

     call put_miller_geometry(tglf_rmin_loc_in, &
          tglf_rmaj_loc_in, &
          tglf_zmaj_loc_in, &
          tglf_drmindx_loc_in, &
          tglf_drmajdx_loc_in, &
          tglf_dzmajdx_loc_in, &
          tglf_kappa_loc_in, &
          tglf_s_kappa_loc_in, &
          tglf_delta_loc_in, &
          tglf_s_delta_loc_in, &
          tglf_zeta_loc_in, &
          tglf_s_zeta_loc_in, &
          tglf_q_loc_in, &
          tglf_q_prime_loc_in, &
          tglf_p_prime_loc_in, &
          tglf_kx0_loc_in)

  elseif (tglf_geometry_flag_in == 2)then

     call put_fourier_geometry(tglf_q_fourier_in,  &
          tglf_q_prime_fourier_in, &
          tglf_p_prime_fourier_in, &
          tglf_nfourier_in, &
          tglf_fourier_in)

  else

     call put_s_alpha_geometry(tglf_rmin_sa_in, &
          tglf_rmaj_sa_in, &
          tglf_q_sa_in, &
          tglf_shat_sa_in, &
          tglf_alpha_sa_in, &
          tglf_xwell_sa_in, &
          tglf_theta0_sa_in, &
          tglf_b_model_sa_in, &
          tglf_ft_model_sa_in)

  endif

  ! Create parameter dump file
  if (tglf_dump_flag_in .eqv. .true.) then
     call tglf_dump_local
     if (tglf_test_flag_in == 1) return
  endif

  if (tglf_use_transport_model_in) then

     !call tglf_tm_mpi
     
     open (unit=14, file="input.dat", action="write")
     
     write (14,*) '1' 

     write (14,"(15(f6.3,x))") tglf_as_in(2), tglf_as_in(3), tglf_betae_in, &
                             tglf_delta_loc_in, tglf_kappa_loc_in, tglf_q_loc_in, &
                             tglf_q_prime_loc_in, tglf_rlns_in(1), tglf_rlts_in(1), &
		             tglf_rlts_in(2), tglf_rmaj_loc_in, tglf_rmin_loc_in, &
		             tglf_s_kappa_loc_in, tglf_taus_in(2), tglf_xnue_in
     
     close(14)
     
     write (*,*) 'RMIN_LOC=', tglf_rmin_loc_in
     write (*,*) 'RMAJ_LOC=', tglf_rmaj_loc_in
     write (*,*) 'KAPPA_LOC=', tglf_kappa_loc_in
     write (*,*) 'S_KAPPA_LOC=', tglf_s_kappa_loc_in
     write (*,*) 'DELTA_LOC=', tglf_delta_loc_in
     write (*,*) 'Q_LOC=', tglf_q_loc_in
     write (*,*) 'Q_PRIME_LOC=', tglf_q_prime_loc_in
     write (*,*) 'XNUE=', tglf_xnue_in
     write (*,*) 'BETAE=', tglf_betae_in
     write (*,*) 'TAUS_2=', tglf_taus_in(2)
     write (*,*) 'AS_2=', tglf_as_in(2)
     write (*,*) 'AS_3=', tglf_as_in(3)
     write (*,*) 'RLNS_1=', tglf_rlns_in(1)
     write (*,*) 'RLTS_1=', tglf_rlts_in(1)
     write (*,*) 'RLTS_2=', tglf_rlts_in(2)
     write (*,*) 'ZEFF=', tglf_zeff_in
     write (*,*) 'S_ZETA_LOC=', tglf_s_zeta_loc_in
     write (*,*) 'S_DELTA_LOC=', tglf_s_delta_loc_in
     write (*,*) 'P_PRIME_LOC=', tglf_p_prime_loc_in
     write (*,*) 'ZMAJ_LOC=', tglf_zmaj_loc_in
     write (*,*) 'DZMAJDX_LOC=', tglf_dzmajdx_loc_in
     write (*,*) 'ZETA_LOC=', tglf_zeta_loc_in
     write (*,*) 'DRMAJDX_LOC=', tglf_drmajdx_loc_in
     
     
     CALL SYSTEM('./run.exe brainfuse_0.net input.dat')
     
     open (unit=15, file="output.avg", action="read")
     
     read(15,*) n, OUT_ENERGY_FLUX_1, OUT_ENERGY_FLUX_3, &
                   OUT_PARTICLE_FLUX_1, OUT_PARTICLE_FLUX_3, &
	           OUT_STRESS_TOR_1, OUT_STRESS_TOR_3   
		   
     write(*,*) 'ELECTRON ENERGY FLUX=', OUT_ENERGY_FLUX_1
     write(*,*) 'ION ENERGY FLUX=', OUT_ENERGY_FLUX_3
     write(*,*) 'ELECTRON PARTICLE FLUX=', OUT_PARTICLE_FLUX_1
     write(*,*) 'ION PARTICLE FLUX=', OUT_PARTICLE_FLUX_3
     write(*,*) 'ELECTRON MOMENTUM FLUX=', OUT_STRESS_TOR_1
     write(*,*) 'ION MOMENTUM FLUX=', OUT_STRESS_TOR_3
               
     close(15)
     
     !---------------------------------------------
     ! Output (normalized to Q_GB)
     ! 
     ! Electrons

     ! Gammae/Gamma_GB
     tglf_elec_pflux_out = OUT_PARTICLE_FLUX_1

     ! Qe/Q_GB
     tglf_elec_eflux_low_out = OUT_ENERGY_FLUX_1
        ! _low already included EM terms
     tglf_elec_eflux_out     = OUT_ENERGY_FLUX_1

     ! Pi_e/Pi_GB
     tglf_elec_mflux_out = OUT_STRESS_TOR_1

     ! S_e/S_GB
     tglf_elec_expwd_out = get_exchange(1,1)  &
          + get_exchange(1,2)                 &
          + get_exchange(1,3)

     ! Ions

     do i_ion=1,5

        ! Gammai/Gamma_GB
        tglf_ion_pflux_out(i_ion) = OUT_PARTICLE_FLUX_3

        ! Qi/Q_GB
        tglf_ion_eflux_low_out(i_ion) = OUT_ENERGY_FLUX_3
           ! _low already included EM terms
        tglf_ion_eflux_out(i_ion)     = OUT_ENERGY_FLUX_3

        ! Pi_i/Pi_GB
        tglf_ion_mflux_out(i_ion) = OUT_STRESS_TOR_3

        ! S_i/S_GB
        tglf_ion_expwd_out(i_ion) = get_exchange(i_ion+1,1)    &
             + get_exchange(i_ion+1,2)                         &
             + get_exchange(i_ion+1,3)

     enddo
     
  else

     ! Run single-ky linear stability
     call tglf_ky

     ! Collect linear eigenvalues
     do n=1,tglf_nmodes_in
        tglf_eigenvalue_out(n) = get_frequency(n) + xi*get_growthrate(n)
     enddo

     ! Print eigenfunction if flag set
     if (iProcTglf == iProc0Tglf .and. tglf_write_wavefunction_flag_in == 1) then
        call write_wavefunction_out(trim(tglf_path_in)//'out.tglf.wavefunction')
     endif

  endif

  ! Diagnostic output
  interchange_DR = get_DR()
  interchange_DM = get_DM()

  call get_error_status(tglf_error_message,tglf_error_status)

end subroutine tglf_run_nn
