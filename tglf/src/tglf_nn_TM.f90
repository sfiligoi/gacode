!-----------------------------------------------------------------
!
     SUBROUTINE tglf_nn_TM
!
!  Execution of the NN
!
     use tglf_pkg
     use tglf_interface
     USE tglf_dimensions
     USE tglf_global
!
     IMPLICIT NONE
!
     INTEGER :: j, is, k
     REAl :: v_bar0, phi_bar0
     REAL :: pflux0(nsm,3),eflux0(nsm,3)
     REAL :: stress_par0(nsm,3),stress_tor0(nsm,3)
     REAL :: exch0(nsm,3)
     REAL :: nsum0(nsm),tsum0(nsm)
!
     REAL, DIMENSION(10) :: tmp
     INTEGER, DIMENSION(10) :: ions_order
!
     real :: OUT_ENERGY_FLUX_1_STD, OUT_ENERGY_FLUX_3_STD, OUT_PARTICLE_FLUX_1_STD
     real :: OUT_PARTICLE_FLUX_3_STD, OUT_STRESS_TOR_3_STD
     real :: OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG
     real :: OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_3_RNG
     integer :: n, j, i
!
! initialize fluxes
!
     do is=1,ns
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
       q_low_out(is) = 0.0
     enddo
!
! '#---------------------------------------------------'
! '# Sort ions by A, Z, Te/Ti'
! '#---------------------------------------------------'
     tmp=0.0
     tmp(1)=1E10
     DO i = 3,tglf_ns_in
         tmp(i)=tglf_mass_in(i)*100000+tglf_zs_in(i)*1000+1./(1+tglf_rlts_in(i))
     ENDDO
     DO i = 1,tglf_ns_in
         k=MAXLOC(tmp, DIM=1)
         ions_order(i)=k
         tmp(k)=0
     ENDDO
!
! Write input file for the NN
!
     open (unit=14, file=TRIM(tglf_path_in)//"input.dat", action="write")
     write (14,*) '1' 
!--------     write (14,"(15(f6.3,x))") tglf_as_in(2), tglf_as_in(ions_order(3)), tglf_betae_in, &
     write (14,"(15(f6.3,x))") tglf_as_in(3), tglf_as_in(2), tglf_betae_in, &
          tglf_delta_loc_in, tglf_kappa_loc_in, tglf_q_loc_in, &
          tglf_q_prime_loc_in, tglf_rlns_in(1), tglf_rlts_in(1), &
          tglf_rlts_in(3), tglf_rmaj_loc_in, tglf_rmin_loc_in, &
          tglf_s_kappa_loc_in, tglf_taus_in(3), tglf_xnue_in
     close(14)
!
! Execute the NN
!
     if (tglf_path_in == "") then
        ! FOR TGLF in standalone mode: no sub-directories
        CALL SYSTEM('/u/ludat/brainfuse_orso/run.exe /u/ludat/testnets/brainfuse_* input.dat')
     else
        ! FOR TGYRO
        CALL SYSTEM('cd '//TRIM(tglf_path_in)//' ;/u/ludat/brainfuse_orso/run.exe /u/ludat/testnets/brainfuse_* input.dat')
     endif
!
! Read outputs
!
     open (unit=15, file=TRIM(tglf_path_in)//"output.avg", action="read")
     read(15,*) n, energy_flux_out(1,1), energy_flux_out(3,1), &
          particle_flux_out(1,1), particle_flux_out(3,1), &
          stress_tor_out(3,1)
     close(15)
!
! Set negative energy fluxes to zero
!    
     if (energy_flux_out(1,1)<0) then
        energy_flux_out(1,1)=0
     endif
     if (energy_flux_out(3,1)<0) then
        energy_flux_out(3,1)=0
     endif
!
! Set Carbon fluxes equal to Deuterium fluxes
!
   !  energy_flux_out(2,1)=energy_flux_out(3,1)
   !  particle_flux_out(2,1)=particle_flux_out(3,1)
   !  stress_tor_out(2,1)=stress_tor_out(3,1)
!
! Read values for checking accuracy of the NN
!
     open (unit=16, file=TRIM(tglf_path_in)//"output.std", action="read")
     open (unit=20, file=TRIM(tglf_path_in)//"output.rng", action="read")
     read(16,*) n, OUT_ENERGY_FLUX_1_STD, OUT_ENERGY_FLUX_3_STD, OUT_PARTICLE_FLUX_1_STD, &
          OUT_PARTICLE_FLUX_3_STD, OUT_STRESS_TOR_3_STD     
     read(20,*) n, OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_3_RNG, OUT_PARTICLE_FLUX_1_RNG, &
          OUT_PARTICLE_FLUX_3_RNG, OUT_STRESS_TOR_3_RNG
     close(16)
     close(20)
!
! Switch between TGLF and the NN depending on 'nn_thrsh' values
!
     if ((OUT_ENERGY_FLUX_1_RNG<tglf_nn_thrsh_energy_in) .and. (OUT_ENERGY_FLUX_3_RNG<tglf_nn_thrsh_energy_in) .and. &
          (OUT_PARTICLE_FLUX_1_RNG<tglf_nn_thrsh_energy_in) .and. (OUT_PARTICLE_FLUX_3_RNG<tglf_nn_thrsh_energy_in) .and. &
          (OUT_STRESS_TOR_3_RNG<tglf_nn_thrsh_energy_in)) then
!
        valid_nn = .TRUE.
!  
        write(*,*) '------------>  THE NN IS RUNNING  <--------------'
!
     else 
!        
        write(*,*) '------------>   TGLF IS RUNNING   <--------------'
!     
     endif
!
      END SUBROUTINE tglf_nn_TM
