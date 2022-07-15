!---------------------------------------------------------
! gftm_run_mpi.f90
!
! PURPOSE:
!  Manage call to gftm simulation for both standalone and
!  TGYRO usage.
!---------------------------------------------------------

#ifdef MPI_GFTM
subroutine gftm_run_mpi()
#else
subroutine gftm_run()
#endif

  use gftm_pkg
  use gftm_interface
  use gftm_global

#ifdef MPI_GFTM
  use gftm_mpi
#endif

  implicit none

  integer :: i_ion,n

  call put_units(gftm_units_in)

  call put_signs(gftm_sign_Bt_in,gftm_sign_It_in)

  call put_species(gftm_ns_in, &
       gftm_zs_in, &
       gftm_mass_in)

  call put_kys(gftm_ky_in)

  call put_gaussian_width(gftm_width_in, &
       gftm_width_min_in, &
       gftm_nwidth_in, &
       gftm_find_width_in)

  call put_eikonal(gftm_new_eikonal_in)

  call put_gradients(gftm_rlns_in, &
       gftm_rlts_in, &
       gftm_vpar_shear_in, &
       gftm_vexb_shear_in)

  call put_averages(gftm_taus_in, &
       gftm_as_in, &
       gftm_vpar_in, &
       gftm_vexb_in, &
       gftm_betae_in, &
       gftm_xnue_in, &
       gftm_zeff_in, &
       gftm_debye_in)

  call put_switches(gftm_iflux_in, &
       gftm_use_bper_in, &
       gftm_use_bpar_in, &
       gftm_use_mhd_rule_in, &
       gftm_use_bisection_in, &
       gftm_use_inboard_detrapped_in, &
       gftm_ibranch_in, &
       gftm_nmodes_in, &
       gftm_nbasis_max_in, &
       gftm_nbasis_min_in, &
       gftm_nxgrid_in, &
       gftm_nky_in, &
       gftm_use_ave_ion_grid_in)

  call put_rare_switches(gftm_theta_trapped_in, &
       gftm_wdia_trapped_in, &
       gftm_park_in, &
       gftm_ghat_in, &
       gftm_gchat_in, &
       gftm_wd_zero_in, &
       gftm_linsker_factor_in, &
       gftm_gradB_factor_in, &
       gftm_filter_in, &
       gftm_damp_psi_in, &
       gftm_damp_sig_in)

  call put_model_parameters(gftm_adiabatic_elec_in, &
       gftm_alpha_e_in, &
       gftm_alpha_p_in, &
       gftm_alpha_mach_in, &
       gftm_alpha_quench_in, &
       gftm_alpha_zf_in, &
       gftm_xnu_factor_in, &
       gftm_debye_factor_in, &
       gftm_etg_factor_in, &
       gftm_rlnp_cutoff_in, &
       gftm_sat_rule_in, &
       gftm_kygrid_model_in, &
       gftm_xnu_model_in, &
       gftm_vpar_model_in, &
       gftm_vpar_shear_model_in)

  if (gftm_geometry_flag_in == 1 ) then

     call put_miller_geometry(gftm_rmin_loc_in, &
          gftm_rmaj_loc_in, &
          gftm_zmaj_loc_in, &
          gftm_drmindx_loc_in, &
          gftm_drmajdx_loc_in, &
          gftm_dzmajdx_loc_in, &
          gftm_kappa_loc_in, &
          gftm_s_kappa_loc_in, &
          gftm_delta_loc_in, &
          gftm_s_delta_loc_in, &
          gftm_zeta_loc_in, &
          gftm_s_zeta_loc_in, &
          gftm_q_loc_in, &
          gftm_q_prime_loc_in, &
          gftm_p_prime_loc_in, &
          gftm_kx0_loc_in)

  elseif (gftm_geometry_flag_in == 2)then

     call put_fourier_geometry(gftm_q_fourier_in,  &
          gftm_q_prime_fourier_in, &
          gftm_p_prime_fourier_in, &
          gftm_nfourier_in, &
          gftm_fourier_in)

  else
     write(*,*)"geometry_flag = ", gftm_geometry_flag_in
     call put_s_alpha_geometry(gftm_rmin_sa_in, &
          gftm_rmaj_sa_in, &
          gftm_q_sa_in, &
          gftm_shat_sa_in, &
          gftm_alpha_sa_in, &
          gftm_xwell_sa_in, &
          gftm_theta0_sa_in, &
          gftm_b_model_sa_in, &
          gftm_ft_model_sa_in)

  endif

  ! Create parameter dump file
  if (gftm_dump_flag_in .eqv. .true.) then
     call gftm_dump_local
     call gftm_dump_global
     if (gftm_test_flag_in == 1) return
  endif

  if (gftm_use_transport_model_in) then

#ifdef MPI_GFTM
        call gftm_tm_mpi
#else
        call gftm_tm
#endif

     !---------------------------------------------
     ! Output (normalized to Q_GB)
     !
     ! Electrons

     ! Gammae/Gamma_GB
     gftm_elec_pflux_out = get_particle_flux(1)

     ! Qe/Q_GB
     gftm_elec_eflux_low_out = get_q_low(1)
     gftm_elec_eflux_out     = get_energy_flux(1)

     ! Pi_e/Pi_GB
     gftm_elec_mflux_out = get_stress_tor(1)

     ! S_e/S_GB
     gftm_elec_expwd_out = get_exchange(1)

     ! Ions

     do i_ion=1,ns_in

             ! Gammai/Gamma_GB
             gftm_ion_pflux_out(i_ion) = get_particle_flux(i_ion+1)
 
             ! Qi/Q_GB
             gftm_ion_eflux_low_out(i_ion) = get_q_low(i_ion+1)
             ! _low already included EM terms
             gftm_ion_eflux_out(i_ion)     = get_energy_flux(i_ion+1)

             ! Pi_i/Pi_GB
             gftm_ion_mflux_out(i_ion) = get_stress_tor(i_ion+1)

             ! S_i/S_GB
             gftm_ion_expwd_out(i_ion) = get_exchange(i_ion+1)

     enddo

     call get_error_status(gftm_error_message,gftm_error_status)

  else

     ! Run single-ky linear stability
     call gftm_ky

     ! Collect linear eigenvalues
     do n=1,gftm_nmodes_in
        gftm_eigenvalue_out(n) = get_frequency(n) + xi*get_growthrate(n)
     enddo

     ! Print eigenfunction if flag set
#ifdef MPI_GFTM
     if (iProcgftm == iProc0gftm .and. gftm_write_wavefunction_flag_in == 1) then
        call write_wavefunction_out(trim(gftm_path_in)//'out.gftm.wavefunction')
     endif
#else
     if (gftm_write_wavefunction_flag_in == 1) then
        call write_wavefunction_out(trim(gftm_path_in)//'out.gftm.wavefunction')
     endif
#endif

  endif

  ! Diagnostic output
  interchange_DR = get_DR()
  interchange_DM = get_DM()

#ifdef MPI_GFTM
end subroutine gftm_run_mpi
#else
end subroutine gftm_run
#endif
