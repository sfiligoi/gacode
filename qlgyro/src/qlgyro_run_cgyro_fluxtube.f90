subroutine qlgyro_run_cgyro_fluxtube

  use mpi
  use qlgyro_globals
  use qlgyro_cgyro_interface
  use tglf_interface
  
  implicit none

  integer :: i_ky_local, i_run, i_colour, i_theta, i_px0_local, i_kypx0
  integer :: ios

  character(len=5) :: kystr, px0str
  real :: ky, timestep, maxtime, kx_max_default, px0

  integer, dimension(:), allocatable :: ky_run, px0_run
  integer, dimension(:), allocatable :: ky_px0_run
  integer :: run_status, restart

  logical :: loop_cycle

  ! Read input.gacode in necessary
  if (cgyro_profile_model_in .eq. 2) then
     ! File needed for output 
     open(unit=1,file=trim(path)//'out.cgyro.info',status='replace')
     call cgyro_make_profiles
     call qlgyro_cgyro_deallocate_arrays
     call map_global2interface
  end if

  ! Maps (Physics) CGYRO interface to QLGYRO globals
  call qlgyro_cgyro_map

  ! Maps (Numerics) CGYRO globals (from input.cgyro file) to CGYRO interface
  call map_global2interface

  ! Maps (Physics) QLGYRO globals to back CGYRO interface
  call cgyro_qlgyro_map

  ! Map all to CGYRO interface to CGYRO globals
  call map_interface2global

  ! Maps QL-GYRO variable to TGLF
  call qlgyro_tglf_map

  call map_global2interface

  call allocate_cgyro_interface
  if (i_proc_global .eq. 0)  call tglf_dump_local

  tglf_nky_in = cgyro_n_toroidal_in - 1

  if (.not. allocated(tglf_ky_spectrum_out)) then
     allocate(tglf_ky_spectrum_out(tglf_nky_in), tglf_dky_spectrum_out(tglf_nky_in))
  end if
  
  do i_ky_local=1, cgyro_n_toroidal_in - 1
     tglf_ky_spectrum_out(i_ky_local) = cgyro_ky_in * i_ky_local
  end do
  
  tglf_dky_spectrum_out = cgyro_ky_in

  call qlgyro_px0_spectrum

  restart = cgyro_restart_flag_in
  tglf_nmodes_in = n_px0

  n_ky = tglf_nky_in

  ! Set up ky px0 grid
  n_kypx0 = n_ky * n_px0
  if (.not. allocated(i_ky)) allocate(i_ky(n_kypx0))
  if (.not. allocated(i_px0)) allocate(i_px0(n_kypx0))

  do i_kypx0=1, n_kypx0
     i_ky(i_kypx0) = (i_kypx0-1) / n_px0 + 1
     i_px0(i_kypx0) = modulo((i_kypx0-1), n_px0) + 1
  end do

  n_thetab = cgyro_n_theta_in * cgyro_n_radial_in
  n_field = cgyro_n_field_in

  ! Eigenvalues
  if (.not.allocated(tglf_eigenvalue_spectrum_out)) allocate(tglf_eigenvalue_spectrum_out(2, tglf_nky_in, tglf_nmodes_in))
  if (.not.allocated(qlgyro_eigenvalue_spectrum_out)) allocate(qlgyro_eigenvalue_spectrum_out(2, n_ky, n_px0))

  ! Fields
  if (.not.allocated(tglf_field_spectrum_out)) allocate(tglf_field_spectrum_out(4, tglf_nky_in, tglf_nmodes_in))
  if (.not.allocated(qlgyro_field_spectrum_out)) allocate(qlgyro_field_spectrum_out(3, n_thetab, n_ky, n_px0))
  if (.not.allocated(qlgyro_theta_ballooning)) allocate(qlgyro_theta_ballooning(n_thetab))
  if (.not.allocated(qlgyro_k_perp)) allocate(qlgyro_k_perp(n_thetab, n_ky, n_px0))
  if (.not.allocated(qlgyro_jacobian)) allocate(qlgyro_jacobian(n_thetab, n_ky, n_px0))

  ! Fluxes
  if (.not.allocated(tglf_flux_spectrum_out)) allocate(tglf_flux_spectrum_out(5, n_species, 3, tglf_nky_in, tglf_nmodes_in))
  if (.not.allocated(qlgyro_flux_spectrum_out)) allocate(qlgyro_flux_spectrum_out(5, n_species, 3, n_ky, n_px0))

  ! Array to check if simulation started for each ky
  if (.not.allocated(ky_run)) allocate(ky_run(n_kypx0))
  if (.not.allocated(ky_color)) allocate(ky_color(n_kypx0))

  ! Grad r geometry factor
  if (.not.allocated(sat_geo_spectrum)) allocate(sat_geo_spectrum(tglf_nky_in))
  if (.not.allocated(kxrms_spectrum)) allocate(kxrms_spectrum(tglf_nky_in))

  ! Initial status of each sim
  ky_run = 0

  ! Initial proc no. for ky run
  ky_color = -1

  ! Initialise status file with all zeros if doesnt exist
  call init_qlgyro_status

  call MPI_BCAST(ky_run, tglf_nky_in, MPI_INTEGER, 0, QLGYRO_COMM_WORLD, ierr)

  tglf_elec_pflux_out = 0.0
  tglf_elec_eflux_out = 0.0
  tglf_elec_mflux_out = 0.0
  tglf_elec_expwd_out = 0.0

  tglf_ion_pflux_out = 0.0
  tglf_ion_eflux_out = 0.0
  tglf_ion_mflux_out = 0.0
  tglf_ion_expwd_out = 0.0

  tglf_particle_flux_out = 0.0
  tglf_energy_flux_out = 0.0
  tglf_stress_tor_out = 0.0
  tglf_exchange_out = 0.0
  tglf_eigenvalue_spectrum_out = 0.0
  tglf_flux_spectrum_out = 0.0
  tglf_field_spectrum_out = 0.0
  qlgyro_eigenvalue_spectrum_out = 0.0
  qlgyro_flux_spectrum_out = 0.0
  qlgyro_field_spectrum_out = 0.0
  qlgyro_k_perp = 0.0
  qlgyro_jacobian = 0.0

  call MPI_barrier(QLGYRO_COMM_WORLD, ierr)

  do i_px0_local=1,n_px0

     ! Initial assignment of jobs to cores
     if (i_px0_local .le. n_parallel .and. i_px0_local .ne. color+1) then
        cycle
     end if

     do i_ky_local=1,tglf_nky_in
        i_kypx0  = (tglf_nky_in * (i_px0_local-1)) + i_ky_local
        call get_qlgyro_status(i_kypx0, run_status, loop_cycle)
     end do
     
     ! If run is being done then skip to next ky
     if (loop_cycle) then
        ! Give some time to allow for other statuses to be set
        call sleep(2)
        call MPI_barrier(CGYRO_COMM_WORLD, ierr)
        cycle
     end if
     run_status = 1

     call write_qlgyro_status(ky_run, ky_color, i_px0_local, run_status, color)

     px0 = px0_spectrum(i_px0_local)

     ! Creates directory for each ky and changes into it for that run
     write(px0str, '(F5.2)') px0
     runpath = trim(trim(path)//"_PX0_"//trim(adjustl(px0str))//"/")

     if (adjoint .eq. 0) then
        call system("mkdir -p "//runpath)
        write(*,*) ' '
        write(*,21) 'Group ',  color, ' running in folder ', trim(runpath), ky_color(i_px0)
     end if

     call MPI_barrier(CGYRO_COMM_WORLD, ierr)

     ! Ensure linear settings for CGYRO
     cgyro_box_size_in = 1     
     cgyro_nonlinear_flag_in = 0
     
     ! Change values after init
     cgyro_px0_in = px0

     cgyro_integration_error_in = 0.0
     cgyro_printout_in = .false.

     cgyro_path_in = trim(runpath)//'/'
     
     call MPI_barrier(CGYRO_COMM_WORLD, ierr)

     cgyro_restart_flag_in = restart

     ! Reset signal to 0
     cgyro_signal_out = 0

     call map_interface2global

     do i_run=1, n_runs
        call cgyro_init_kernel
        call cgyro_kernel
        call qlgyro_cgyro_cleanup
        call cgyro_final_kernel
        run_status = 2
     end do

     do i_ky_local=1, tglf_nky_in

        ky = tglf_ky_spectrum_out(i_ky_local)

        !if (abs(cgyro_omega_error_out(i_ky_local+1)) .lt. cgyro_freq_tol_in) then
        !   run_status = 3
        !end if

        ! EAB: note about cgyro_omega_error_out
        ! cgyro_omega_error_out and freq_err were removed from CGYRO
        ! cgyro_omega_error_out is not computed in cgyro for each toroidal mode
        ! (it is a function of n_toroidal_per_proc, not n_toroidal
        ! so would only be computed for each ky if n_toroidal_per_proc = n_toroidal
        ! in the parallel layout at runtime)
        ! there needs to be a different way to check convergence
        ! it would be more efficient to just run the cgyro's at each ky separately
        ! to convergence (e.g. via tgyro)
        ! for now, we will just assume the code is converged at each ky
        run_status = 3
        
        i_kypx0  = (tglf_nky_in * i_px0_local) + i_ky_local

        call write_qlgyro_status(ky_run, ky_color, i_kypx0, run_status, color)

        if (run_status .eq. 3) then

           if (adjoint == 0 ) write(*,22) 'CGYRO run converged at ky/px0 =', ky,'/',px0

           ! CGYRO eigenvalues
           qlgyro_eigenvalue_spectrum_out(2,:,i_px0_local) = real(cgyro_omega_out(i_ky_local))
           qlgyro_eigenvalue_spectrum_out(1,:,i_px0_local) = aimag(cgyro_omega_out(i_ky_local))
           
           ! CGYRO has electrons at last index, TGLF has them in first so load them seperately
           qlgyro_flux_spectrum_out(1:3, 1, :cgyro_n_field_in, i_ky_local, i_px0_local) =         &
                cgyro_gbflux_out(1:3, elec_ind, :cgyro_n_field_in, i_ky_local) * ky * i_ky_local * (-cgyro_btccw_in)
           
           qlgyro_flux_spectrum_out(1:3, 2:, :cgyro_n_field_in, i_ky_local, i_px0_local) = &
                cgyro_gbflux_out(1:3, ion_start:ion_end, :cgyro_n_field_in, i_ky_local) * ky * i_ky_local * (-cgyro_btccw_in)
           
           qlgyro_theta_ballooning = cgyro_thetab_out
           qlgyro_field_spectrum_out(:, :, i_ky_local, i_px0_local) = cgyro_wavefunction_out(i_ky_local, :, :)
           
           qlgyro_k_perp(:, i_ky_local, i_px0_local) = cgyro_k_perp_out(i_ky_local, :)
           qlgyro_jacobian(:, i_ky_local, i_px0_local) = cgyro_jacobian_out
           
        ! If run hasn't converged, sent eigenvalue/growth to 0
        else if (run_status .eq. 2) then
           if (adjoint .eq. 0) write(*,22) 'CGYRO not converged at ky/px0 = ', ky,'/',px0
           qlgyro_eigenvalue_spectrum_out(:,i_ky_local,i_px0_local) = 0.0
           qlgyro_flux_spectrum_out(1:4, :, :cgyro_n_field_in, i_ky_local, i_px0_local) = 0.0
        end if

     end do
     call MPI_BARRIER(CGYRO_COMM_WORLD, ierr)
  end do

  call MPI_BARRIER(QLGYRO_COMM_WORLD, ierr)
  
  call read_qlgyro_status(ky_run, ky_color)
 
  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1, *) 'INFO: Ran all kys'
     write(1, *) 'B_unit = ', bunit
     close(1)
  endif

  
21 format(A6, I3, A24, A19, I3)
22 format(A32, F5.2, A1, F5.2)
end subroutine qlgyro_run_cgyro_fluxtube

  
    
