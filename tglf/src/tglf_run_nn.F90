!---------------------------------------------------------
! tglf_run_mpi.f90
!
! PURPOSE:
!  Manage call to TGLF simulation for both standalone and
!  TGYRO usage.
!---------------------------------------------------------

#ifdef MPI_TGLF
subroutine tglf_run_nn_mpi()
#else
subroutine tglf_run_nn()
#endif

  use tglf_pkg
  use tglf_interface
  
  #ifdef MPI_TGLF
  use tglf_mpi
  #endif
  
  implicit none

  real :: OUT_ENERGY_FLUX_1, OUT_PARTICLE_FLUX_1, OUT_STRESS_TOR_1
  real :: OUT_ENERGY_FLUX_3, OUT_PARTICLE_FLUX_3, OUT_STRESS_TOR_3
  real :: OUT_ENERGY_FLUX_1_STD, OUT_ENERGY_FLUX_3_STD, OUT_PARTICLE_FLUX_1_STD
  real :: OUT_PARTICLE_FLUX_3_STD, OUT_STRESS_TOR_1_STD, OUT_STRESS_TOR_3_STD
  real :: AS_2_LIM, AS_3_LIM, BETAE_LIM, DELTA_LOC_LIM, KAPPA_LOC_LIM, Q_LOC_LIM
  real :: Q_PRIME_LOC_LIM, RLNS_1_LIM, RLTS_1_LIM, RLTS_2_LIM, RMAJ_LOC_LIM 
  real :: RMIN_LOC_LIM, S_KAPPA_LOC_LIM, TAUS_2_LIM, XNUE_LIM
  real :: OUT_ENERGY_FLUX_1_LIM, OUT_ENERGY_FLUX_3_LIM, OUT_PARTICLE_FLUX_1_LIM
  real :: OUT_PARTICLE_FLUX_3_LIM, OUT_STRESS_TOR_1_LIM, OUT_STRESS_TOR_3_LIM
  real :: OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG
  real :: OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_1_RNG, OUT_STRESS_TOR_3_RNG

  real :: nn_in_lim=10000
  integer :: nn_check=0
    
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

     #ifdef MPI_TGLF
     ! Check proper processor setup
     if (iProcTglf .ne. iProc0Tglf) then
        write (*,*) 'wrong processor setup: parallelization not yet supported for NN'
        stop
     endif
     #endif

     ! Write input file for the NN
     open (unit=14, file=TRIM(tglf_path_in)//"input.dat", action="write")
     write (14,*) '1' 
     write (14,"(15(f6.3,x))") tglf_as_in(3), tglf_as_in(2), tglf_betae_in, &
          tglf_delta_loc_in, tglf_kappa_loc_in, tglf_q_loc_in, &
          tglf_q_prime_loc_in, tglf_rlns_in(1), tglf_rlts_in(1), &
          tglf_rlts_in(3), tglf_rmaj_loc_in, tglf_rmin_loc_in, &
          tglf_s_kappa_loc_in, tglf_taus_in(3), tglf_xnue_in
     close(14)

     ! Execute the NN
     if (tglf_path_in == "") then
        ! FOR TGLF    O N L Y
        write(*,*) 'tglf_run_nn --> tglf_path_in is blank'
        CALL SYSTEM('/u/ludat/brainfuse_orso/run.exe /u/ludat/testnets/brainfuse_* input.dat')
     else
        write(*,*) 'tglf_run_nn --> tglf_path_in: ', trim(tglf_path_in)
        ! FOR TGYRO
        CALL SYSTEM('cd '//TRIM(tglf_path_in)//' ;/u/ludat/brainfuse_orso/run.exe /u/ludat/testnets/brainfuse_* input.dat')
     endif

     ! Read outputs
     open (unit=15, file=TRIM(tglf_path_in)//"output.avg", action="read")
     read(15,*) n, OUT_ENERGY_FLUX_1, OUT_ENERGY_FLUX_3, &
          OUT_PARTICLE_FLUX_1, OUT_PARTICLE_FLUX_3, &
          OUT_STRESS_TOR_1, OUT_STRESS_TOR_3
     close(15)

     ! Read values for checking accuracy of the NN
     open (unit=16, file=TRIM(tglf_path_in)//"output.std", action="read")
     open (unit=17, file=TRIM(tglf_path_in)//"input.lim", action="read")
     open (unit=19, file=TRIM(tglf_path_in)//"output.lim", action="read")
     open (unit=20, file=TRIM(tglf_path_in)//"output.rng", action="read")
     read(16,*) n, OUT_ENERGY_FLUX_1_STD, OUT_ENERGY_FLUX_3_STD, OUT_PARTICLE_FLUX_1_STD, &
          OUT_PARTICLE_FLUX_3_STD, OUT_STRESS_TOR_1_STD, OUT_STRESS_TOR_3_STD     
     read(17,*) n, AS_2_LIM, AS_3_LIM, BETAE_LIM, DELTA_LOC_LIM, KAPPA_LOC_LIM, Q_LOC_LIM, &
          Q_PRIME_LOC_LIM, RLNS_1_LIM, RLTS_1_LIM, RLTS_2_LIM, RMAJ_LOC_LIM, &
          RMIN_LOC_LIM, S_KAPPA_LOC_LIM, TAUS_2_LIM, XNUE_LIM
     read(19,*) n, OUT_ENERGY_FLUX_1_LIM, OUT_ENERGY_FLUX_3_LIM, OUT_PARTICLE_FLUX_1_LIM, &
          OUT_PARTICLE_FLUX_3_LIM, OUT_STRESS_TOR_1_LIM, OUT_STRESS_TOR_3_LIM
     read(20,*) n, OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG, &
          OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_1_RNG, OUT_STRESS_TOR_3_RNG
     close(16)
     close(17)
     close(19)
     close(20)

     ! If only 1 NN: check 'input.lim', else: check JUST 'output.rng'
     if ((OUT_ENERGY_FLUX_1_RNG==0) .and. (OUT_ENERGY_FLUX_3_RNG==0) .and. (OUT_PARTICLE_FLUX_1_RNG==0) .and. &
          (OUT_PARTICLE_FLUX_3_RNG==0) .and. (OUT_STRESS_TOR_1_RNG==0) .and. (OUT_STRESS_TOR_3_RNG==0)) then
        nn_in_lim=2
     else 
        nn_in_lim=9999 
     endif

     ! Switch between TGLF and the NN depending on 'input.lim' values
     if ((abs(AS_2_LIM)<nn_in_lim) .and. (abs(AS_3_LIM)<nn_in_lim) .and. (abs(BETAE_LIM)<nn_in_lim) .and. &
          (abs(DELTA_LOC_LIM)<nn_in_lim) .and. (abs(KAPPA_LOC_LIM)<nn_in_lim) .and. (abs(Q_LOC_LIM)<nn_in_lim) .and. &
          (abs(Q_PRIME_LOC_LIM)<nn_in_lim) .and. (abs(RLNS_1_LIM)<nn_in_lim) .and. (abs(RLTS_1_LIM)<nn_in_lim) .and. &
          (abs(RLTS_2_LIM)<nn_in_lim) .and. (abs(RMAJ_LOC_LIM)<nn_in_lim) .and. (abs(RMIN_LOC_LIM)<nn_in_lim) .and. &
          (abs(S_KAPPA_LOC_LIM)<nn_in_lim) .and. (abs(TAUS_2_LIM)<nn_in_lim) .and. (abs(XNUE_LIM)<nn_in_lim) .and. &
          (OUT_ENERGY_FLUX_1_RNG<tglf_nn_thrsh_energy_in) .and. (OUT_ENERGY_FLUX_3_RNG<tglf_nn_thrsh_energy_in) .and. &
          (OUT_PARTICLE_FLUX_1_RNG<tglf_nn_thrsh_energy_in) .and. (OUT_PARTICLE_FLUX_3_RNG<tglf_nn_thrsh_energy_in) .and. &
          (OUT_STRESS_TOR_1_RNG<tglf_nn_thrsh_energy_in) .and. (OUT_STRESS_TOR_3_RNG<tglf_nn_thrsh_energy_in)) then
	  
	write(*,*) '------------>  THE NN IS RUNNING  <--------------'
        write(*,"(6(f6.3,x))") OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG, &
             OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_1_RNG, OUT_STRESS_TOR_3_RNG

	! NN case: store results in the file 'results.nn'
        open (unit=18, file=TRIM(tglf_path_in)//"results.nn", position="append", action="write")
        write (18,*) tglf_path_in      
        write (18,*) 'STD: OEF1  OEF3  OPF1  OPF3  OST1  OST3'
        write (18,"(6(f6.3,x))") OUT_ENERGY_FLUX_1_STD, OUT_ENERGY_FLUX_3_STD, OUT_PARTICLE_FLUX_1_STD, &
             OUT_PARTICLE_FLUX_3_STD, OUT_STRESS_TOR_1_STD, OUT_STRESS_TOR_3_STD
        write (18,*) 'NN_OUT: OEF1  OEF3  OPF1  OPF3  OST1  OST3'
        write (18,"(6(f6.3,x))") OUT_ENERGY_FLUX_1, OUT_ENERGY_FLUX_3, OUT_PARTICLE_FLUX_1, &
             OUT_PARTICLE_FLUX_3, OUT_STRESS_TOR_1, OUT_STRESS_TOR_3
        write (18,*) 'LIM: AS2  AS3  BETAE  DELTA  KAPPA  Q  QPRIME  RLNS1  RLTS1  RLTS2  RMAJ  RMIN  SKAPPA  TAUS2  XNUE'
        write (18,"(15(f8.3,x))") AS_2_LIM, AS_3_LIM, BETAE_LIM, DELTA_LOC_LIM, KAPPA_LOC_LIM, Q_LOC_LIM, &
             Q_PRIME_LOC_LIM, RLNS_1_LIM, RLTS_1_LIM, RLTS_2_LIM, RMAJ_LOC_LIM, &
             RMIN_LOC_LIM, S_KAPPA_LOC_LIM, TAUS_2_LIM, XNUE_LIM
        write (18,*) 'LIM: OEF1  OEF3  OPF1  OPF3  OST1  OST3'
        write (18,"(6(f6.3,x))") OUT_ENERGY_FLUX_1_LIM, OUT_ENERGY_FLUX_3_LIM, OUT_PARTICLE_FLUX_1_LIM, &
             OUT_PARTICLE_FLUX_3_LIM, OUT_STRESS_TOR_1_LIM, OUT_STRESS_TOR_3_LIM
        close(18)

        ! Write outputs for all the radial points where the NN is used
        open (unit=25, file="results.nns", position="append", action="write")
        write (25,*) tglf_path_in
        write (25,*) 'NN: OEF1  OEF3  OPF1  OPF3  OST1  OST3'
        write (25,"(6(f6.3,x))") OUT_ENERGY_FLUX_1, OUT_ENERGY_FLUX_3, OUT_PARTICLE_FLUX_1, &
             OUT_PARTICLE_FLUX_3, OUT_STRESS_TOR_1, OUT_STRESS_TOR_3
        write (25,*)
        close(25)

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
        !tglf_elec_expwd_out = get_exchange(1,1)  &
        !     + get_exchange(1,2)                 &
        !     + get_exchange(1,3)

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
           !tglf_ion_expwd_out(i_ion) = get_exchange(i_ion+1,1)    &
           !     + get_exchange(i_ion+1,2)                         &
           !     + get_exchange(i_ion+1,3)

        enddo


     else

        ! TGLF default case
        write(*,*) '------------>  ___ TGLF IS RUNNING ___  <--------------'
        write(*,*) 'r=', tglf_rmin_loc_in, 'r lim=', RMIN_LOC_LIM
        write(*,*) 'R=', tglf_rmaj_loc_in, 'R lim=', RMAJ_LOC_LIM
        write(*,*) 'beta=', tglf_betae_in, 'betae lim=', BETAE_LIM
        write(*,*) 'delta=', tglf_delta_loc_in, 'delta lim=', DELTA_LOC_LIM
        write(*,*) 'kappa=', tglf_kappa_loc_in, 'kappa lim=', KAPPA_LOC_LIM
        write(*,*) 'delta=', tglf_delta_loc_in, 'delta lim=', DELTA_LOC_LIM
        write(*,*) 'q=', tglf_q_loc_in, 'q lim=', Q_LOC_LIM
        write(*,*) 'q prime=', tglf_q_prime_loc_in, 'q prime lim=', Q_PRIME_LOC_LIM
        write(*,*) 's kappa=', tglf_s_kappa_loc_in, 's kapa lim=', S_KAPPA_LOC_LIM
        write(*,*) 'AAAs2=', tglf_as_in(2), 'as2 lim=', AS_2_LIM
        write(*,*) 'as3=', tglf_as_in(3), 'as3 lim=', AS_3_LIM
        write(*,*) 'xnue=', tglf_xnue_in, 'xnue lim=', XNUE_LIM
        write(*,*) 'taus2=', tglf_taus_in(3), 'taus2 lim=', TAUS_2_LIM
        write(*,*) 'rlns1=', tglf_rlns_in(1), 'rlns1 lim=', RLNS_1_LIM
        write(*,*) 'rlts1=', tglf_rlts_in(1), 'rlts1 lim=', RLTS_1_LIM
        write(*,*) 'rlts2=', tglf_rlts_in(3), 'r lim=', RLTS_2_LIM
        write(*,"(6(f6.3,x))") OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG, &
             OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_1_RNG, OUT_STRESS_TOR_3_RNG

        
	#ifdef MPI_TGLF
	call tglf_tm_mpi
	write(*,*) 'hello'
	#else
	call tglf_tm
        write(*,*) 'no hello'
	#endif

        !---------------------------------------------
        ! Output (normalized to Q_GB)
        ! 
        ! Electrons

        ! Gammae/Gamma_GB
        tglf_elec_pflux_out = get_particle_flux(1,1)  &
             + get_particle_flux(1,2)                 &
             + get_particle_flux(1,3)

        write(*,*) 'e particle flux for TLGF and NN: ', tglf_elec_pflux_out, OUT_PARTICLE_FLUX_1

        ! Qe/Q_GB
        tglf_elec_eflux_low_out = get_q_low(1)
        tglf_elec_eflux_out     = get_energy_flux(1,1) &
             + get_energy_flux(1,2)                    &
             + get_energy_flux(1,3)

        write(*,*) 'e energy flux for TLGF and NN: ', tglf_elec_eflux_out, OUT_ENERGY_FLUX_1

        ! Pi_e/Pi_GB
        tglf_elec_mflux_out = get_stress_tor(1,1)   &
             + get_stress_tor(1,2)                  &
             + get_stress_tor(1,3)

        ! S_e/S_GB
        tglf_elec_expwd_out = get_exchange(1,1)  &
             + get_exchange(1,2)                 &
             + get_exchange(1,3)

        ! Ions

        do i_ion=1,5

           ! Gammai/Gamma_GB
           tglf_ion_pflux_out(i_ion) = get_particle_flux(i_ion+1,1) &
                + get_particle_flux(i_ion+1,2)                      &
                + get_particle_flux(i_ion+1,3)

           ! Qi/Q_GB
           tglf_ion_eflux_low_out(i_ion) = get_q_low(i_ion+1)
           ! _low already included EM terms
           tglf_ion_eflux_out(i_ion)     = get_energy_flux(i_ion+1,1) &
                + get_energy_flux(i_ion+1,2)                          &
                + get_energy_flux(i_ion+1,3)

           ! Pi_i/Pi_GB
           tglf_ion_mflux_out(i_ion) = get_stress_tor(i_ion+1,1)  &
                + get_stress_tor(i_ion+1,2)                       &
                + get_stress_tor(i_ion+1,3)

           ! S_i/S_GB
           tglf_ion_expwd_out(i_ion) = get_exchange(i_ion+1,1)    &
                + get_exchange(i_ion+1,2)                         &
                + get_exchange(i_ion+1,3)

        enddo


        ! TGLF case: store results in the file 'results.tglf'
        open (unit=21, file=TRIM(tglf_path_in)//"results.tglf", position="append", action="write")
        write (21,*) tglf_path_in
        write (21,*) 'TGLF: OEF1  OEF3  OPF1  OPF3  OST1  OST3'
        write (21,"(6(f6.3,x))") tglf_elec_eflux_out, tglf_ion_eflux_out(1), tglf_elec_pflux_out, &
             tglf_ion_pflux_out(1), tglf_elec_mflux_out, tglf_ion_mflux_out(1)
        close(21)

        ! Write outputs for all the radial points where TGLF is used
        open (unit=27, file="results.nns", position="append", action="write")
        write (27,*) tglf_path_in
        write (27,*) 'TGLF: OEF1  OEF3  OPF1  OPF3  OST1  OST3'
        write (27,"(6(f6.3,x))") tglf_elec_eflux_out, tglf_ion_eflux_out(1), tglf_elec_pflux_out, &
             tglf_ion_pflux_out(1), tglf_elec_mflux_out, tglf_ion_mflux_out(1)
        write (27,*)
        close(27)


        call get_error_status(tglf_error_message,tglf_error_status)

        IF (tglf_error_status.EQ.0) THEN
           CALL tglf_harvest_local
        ENDIF


     end if

  else

     ! Run single-ky linear stability
     call tglf_ky

     ! Collect linear eigenvalues
     do n=1,tglf_nmodes_in
        tglf_eigenvalue_out(n) = get_frequency(n) + xi*get_growthrate(n)
     enddo

     ! Print eigenfunction if flag set
     #ifdef MPI_TGLF
     if (iProcTglf == iProc0Tglf .and. tglf_write_wavefunction_flag_in == 1) then
        call write_wavefunction_out(trim(tglf_path_in)//'out.tglf.wavefunction')
     endif     
     #else
     if (tglf_write_wavefunction_flag_in == 1) then
        call write_wavefunction_out(trim(tglf_path_in)//'out.tglf.wavefunction')
     endif
     #endif

  endif

  ! Diagnostic output
  interchange_DR = get_DR()
  interchange_DM = get_DM()

#ifdef MPI_TGLF
end subroutine tglf_run_nn_mpi
#else
end subroutine tglf_run_nn
#endif
