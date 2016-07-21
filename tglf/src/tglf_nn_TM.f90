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
     character(len=1000) :: nn_executable
     character(len=1000) :: nn_files

     integer :: j,is,ierr
     real :: v_bar0, phi_bar0
     real :: pflux0(nsm,3),eflux0(nsm,3)
     real :: stress_par0(nsm,3),stress_tor0(nsm,3)
     real :: exch0(nsm,3)
     real :: nsum0(nsm),tsum0(nsm)
!
     real, DIMENSION(10) :: tmp
     integer, DIMENSION(10) :: ions_order
!
     real :: OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_i_RNG
     real :: OUT_PARTICLE_FLUX_1_RNG, OUT_STRESS_TOR_i_RNG
     integer :: n, i

     include 'brainfuse_lib.inc'
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
     tmp(2)=1E9
     DO i = 3,tglf_ns_in
         tmp(i)=tglf_mass_in(i)*100000+tglf_zs_in(i)*1000+1./(1+tglf_rlts_in(i))
     ENDDO
     DO i = 1,tglf_ns_in
         j=MAXLOC(tmp, DIM=1)
         ions_order(i)=j
         tmp(j)=0
     ENDDO

    if ( (tglf_mass_in(2) .ne. 1) .or. &
         (tglf_zs_in(2)   .ne. 1) .or. &
         (tglf_mass_in(3) .ne. 6) .or. &
         (tglf_zs_in(3)   .ne. 6) ) then
        write(*,*)'A',tglf_mass_in
        write(*,*)'Z',tglf_zs_in
        write (*,*)'NN trained only with D C ions'
        stop
    endif

! Write input file for the NN
!
     open (unit=14, file=TRIM(tglf_path_in)//"input.dat", action="write")
     write (14,*) '1'
     write (14,"(15(f6.3,1x))") tglf_as_in(2),       & ! AS_2
                                tglf_as_in(3),       & ! AS_3
                                tglf_betae_in,       & ! BETAE
                                tglf_delta_loc_in,   & ! DELTA_LOC
                                tglf_kappa_loc_in,   & ! KAPPA_LOC
                                tglf_q_loc_in,       & ! Q_LOC
                                tglf_q_prime_loc_in, & ! Q_PRIME_LOC
                                tglf_rlns_in(1),     & ! RLNS_1
                                tglf_rlns_in(2),     & ! RLNS_2
                                tglf_rlns_in(3),     & ! RLNS_3
                                tglf_rlts_in(1),     & ! RLTS_1
                                tglf_rlts_in(2),     & ! RLTS_2
                                tglf_rmaj_loc_in,    & ! RMAJ_LOC
                                tglf_rmin_loc_in,    & ! RMIN_LOC
                                tglf_s_kappa_loc_in, & ! S_KAPPA_LOC
                                tglf_taus_in(2),     & ! TAUS_2
                                tglf_xnue_in           ! XNUE
     close(14)
!
! Execute the NN
!
     call get_environment_variable('TGLFNN_EXEC',nn_executable)
     call get_environment_variable('TGLFNN_MODEL',nn_files)
     if (len(trim(nn_executable))==0) then
        write(*,*)'TGLFNN_EXEC environmental variable must be defined to use NN model'
        stop
     endif
     if (len(trim(trim(nn_files)))==0) then
        write(*,*)'TGLFNN_MODEL environmental variable must be defined to use NN model'
        stop
     endif
     if (tglf_path_in .ne. "") then
        nn_executable='cd '//TRIM(tglf_path_in)//' ;'//trim(nn_executable)
     endif
     call gacode_system(trim(nn_executable)//' '//trim(nn_files)//' input.dat')
!
! Read outputs
!
     open (unit=15, file=TRIM(tglf_path_in)//"output.avg", action="read")
     read(15,*) n,  energy_flux_out(1,1),   &
                    energy_flux_out(3,1),   &
                    particle_flux_out(1,1), &
                    stress_tor_out(3,1)
     close(15)
!
! Read values for checking accuracy of the NN
!
     open (unit=15, file=TRIM(tglf_path_in)//"output.rng", action="read")
     read(15,*) n,  OUT_ENERGY_FLUX_1_RNG,   &
                    OUT_ENERGY_FLUX_i_RNG,   &
                    OUT_PARTICLE_FLUX_1_RNG, &
                    OUT_STRESS_TOR_i_RNG
     close(15)
!
! Switch between TGLF and the NN depending on 'nn_max_error' values
!
     if ((OUT_ENERGY_FLUX_1_RNG   < tglf_nn_max_error_in) .and. &
         (OUT_ENERGY_FLUX_i_RNG   < tglf_nn_max_error_in) .and. &
         (OUT_PARTICLE_FLUX_1_RNG < tglf_nn_max_error_in) .and. &
         (OUT_STRESS_TOR_i_RNG    < tglf_nn_max_error_in)) then
        valid_nn = .TRUE.
        write(*,*) '============>    NN    RUNNING    <=============='
     else
        write(*,*) '------------>   TGLF    RUNNING   <--------------'
     endif
!
      END SUBROUTINE tglf_nn_TM
