!-----------------------------------------------------------------
!
      SUBROUTINE tglf_nn_TM
!
!  Main transport model subroutine.
!  Calls linear TGLF over a spectrum of ky's and computes spectral integrals of 
!  field, intensity and fluxes.
!
      use tglf_pkg
      use tglf_interface

      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum

      IMPLICIT NONE


      INTEGER :: i,j,is,imax
      REAL :: dky
      REAL :: phi_bar0,phi_bar1
      REAl :: v_bar0,v_bar1
      REAL :: dky0,dky1,ky0,ky1
      REAL :: pflux0(nsm,3),eflux0(nsm,3)
      REAL :: stress_par0(nsm,3),stress_tor0(nsm,3)
      REAL :: exch0(nsm,3)
      REAL :: nsum0(nsm),tsum0(nsm)
      REAL :: pflux1(nsm,3),eflux1(nsm,3)
      REAL :: stress_par1(nsm,3),stress_tor1(nsm,3)
      REAL :: exch1(nsm,3)
      REAL :: nsum1(nsm),tsum1(nsm)


  real :: OUT_ENERGY_FLUX_1_STD, OUT_ENERGY_FLUX_3_STD, OUT_PARTICLE_FLUX_1_STD
  real :: OUT_PARTICLE_FLUX_3_STD, OUT_STRESS_TOR_1_STD, OUT_STRESS_TOR_3_STD

  real :: OUT_ENERGY_FLUX_1_LIM, OUT_ENERGY_FLUX_3_LIM, OUT_PARTICLE_FLUX_1_LIM
  real :: OUT_PARTICLE_FLUX_3_LIM, OUT_STRESS_TOR_1_LIM, OUT_STRESS_TOR_3_LIM

  real :: OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG
  real :: OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_1_RNG, OUT_STRESS_TOR_3_RNG

  real :: AS_2_LIM, AS_3_LIM, BETAE_LIM, DELTA_LOC_LIM, KAPPA_LOC_LIM, Q_LOC_LIM
  real :: Q_PRIME_LOC_LIM, RLNS_1_LIM, RLTS_1_LIM, RLTS_2_LIM, RMAJ_LOC_LIM 
  real :: RMIN_LOC_LIM, S_KAPPA_LOC_LIM, TAUS_2_LIM, XNUE_LIM

  integer :: n


  real :: nn_in_lim=10000
  integer :: nn_check=0


!
! initialize fluxes
!
      do is=ns0,ns
        do j=1,3 
          particle_flux_out(is,j) = 0.0
          energy_flux_out(is,j) = 0.0
          stress_par_out(is,j) = 0.0
          stress_tor_out(is,j) = 0.0
          exchange_out(is,j) = 0.0
          pflux0(is,j) = 0.0
          eflux0(is,j) = 0.0
          stress_par0(is,j) = 0.0
          stress_tor0(is,j) = 0.0
          exch0(is,j) = 0.0
        enddo
        n_bar_sum_out(is) = 0.0
        t_bar_sum_out(is) = 0.0
        q_low_out(is) = 0.0
        nsum0(is) = 0.0
        tsum0(is) = 0.0
      enddo
      phi_bar_sum_out = 0.0
      v_bar_sum_out = 0.0
      v_bar0 = 0.0
      phi_bar0 = 0.0     
!

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
        write(*,*) 'tglf_run --> tglf_path_in is blank'
        CALL SYSTEM('/u/ludat/brainfuse_orso/run.exe /u/ludat/testnets/brainfuse_* input.dat')
     else
        write(*,*) 'tglf_run --> tglf_path_in: ', trim(tglf_path_in)
        ! FOR TGYRO
        CALL SYSTEM('cd '//TRIM(tglf_path_in)//' ;/u/ludat/brainfuse_orso/run.exe /u/ludat/testnets/brainfuse_* input.dat')
     endif

     ! Read outputs
     open (unit=15, file=TRIM(tglf_path_in)//"output.avg", action="read")
     read(15,*) n, energy_flux_out(1,1), energy_flux_out(3,1), &
          particle_flux_out(1,1), particle_flux_out(3,1), &
          stress_tor_out(1,1), stress_tor_out(3,1)
     close(15)
     
     energy_flux_out(2,1)=energy_flux_out(3,1)
     particle_flux_out(2,1)=particle_flux_out(3,1)
     stress_tor_out(2,1)=stress_tor_out(3,1)
     

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

        valid_nn = .TRUE.
	  
        write(*,*) '------------>  THE NN IS RUNNING  <--------------'
        write(*,"(6(f6.3,x))") OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG, &
             OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_1_RNG, OUT_STRESS_TOR_3_RNG
	     
     else 
        
        write(*,*) '------------>  ___ TGLF IS RUNNING ___  <--------------'
     
     endif
!
      END SUBROUTINE tglf_nn_TM

        ! TGLF default case
       ! write(*,*) '------------>  ___ TGLF IS RUNNING ___  <--------------'
       ! write(*,*) 'r=', tglf_rmin_loc_in, 'r lim=', RMIN_LOC_LIM
       ! write(*,*) 'R=', tglf_rmaj_loc_in, 'R lim=', RMAJ_LOC_LIM
       ! write(*,*) 'beta=', tglf_betae_in, 'betae lim=', BETAE_LIM
       ! write(*,*) 'delta=', tglf_delta_loc_in, 'delta lim=', DELTA_LOC_LIM
       ! write(*,*) 'kappa=', tglf_kappa_loc_in, 'kappa lim=', KAPPA_LOC_LIM
       ! write(*,*) 'delta=', tglf_delta_loc_in, 'delta lim=', DELTA_LOC_LIM
       ! write(*,*) 'q=', tglf_q_loc_in, 'q lim=', Q_LOC_LIM
       ! write(*,*) 'q prime=', tglf_q_prime_loc_in, 'q prime lim=', Q_PRIME_LOC_LIM
       ! write(*,*) 's kappa=', tglf_s_kappa_loc_in, 's kapa lim=', S_KAPPA_LOC_LIM
       ! write(*,*) 'AAAs2=', tglf_as_in(2), 'as2 lim=', AS_2_LIM
       ! write(*,*) 'as3=', tglf_as_in(3), 'as3 lim=', AS_3_LIM
       ! write(*,*) 'xnue=', tglf_xnue_in, 'xnue lim=', XNUE_LIM
       ! write(*,*) 'taus2=', tglf_taus_in(3), 'taus2 lim=', TAUS_2_LIM
       ! write(*,*) 'rlns1=', tglf_rlns_in(1), 'rlns1 lim=', RLNS_1_LIM
       ! write(*,*) 'rlts1=', tglf_rlts_in(1), 'rlts1 lim=', RLTS_1_LIM
       ! write(*,*) 'rlts2=', tglf_rlts_in(3), 'r lim=', RLTS_2_LIM
       ! write(*,"(6(f6.3,x))") OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG, &
       !      OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_1_RNG, OUT_STRESS_TOR_3_RNG

