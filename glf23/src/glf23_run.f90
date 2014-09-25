!---------------------------------------------------------
! glf23_run.f90
!
! PURPOSE:
!  Manage call to glf23 simulation for both standalone and
!  TGYRO usage.
!---------------------------------------------------------

subroutine glf23_run()

  use glf23_interface

  implicit none

  integer :: i_ion,n
  complex :: xi=(0.0,1.0)


  call put_model_parameters(glf23_alpha_p_in, &
       glf23_alpha_quench_in, &
       glf23_version_in)

  call put_species(glf23_ns_in, &
       glf23_zs_in, &
       glf23_mass_in)

  call put_kys(glf23_ky_in)

  call put_gradients(glf23_rlns_in, &
       glf23_rlts_in, &
       glf23_vpar_shear_in, &
       glf23_vexb_shear_in)

  call put_averages(glf23_taus_in, &
       glf23_as_in, &
       glf23_betae_in, &
       glf23_xnue_in)

  call put_s_alpha_geometry(glf23_rmin_sa_in, &
       glf23_rmaj_sa_in, &
       glf23_q_sa_in, &
       glf23_shat_sa_in, &
       glf23_alpha_sa_in, &
       glf23_xwell_sa_in, &
       glf23_theta0_sa_in)

  ! Create paramter dump files
  if (glf23_dump_flag_in .eqv. .true.) then
     call glf23_dump_local
  endif

  if (glf23_use_transport_model_in) then

     call glf2d

     !---------------------------------------------
     ! Output (normalized to Q_GB)
     ! 
     ! Electrons

     ! Gammae/Gamma_GB
     glf23_elec_pflux_out = get_particle_flux(1,1)

     ! Qe/Q_GB
     glf23_elec_eflux_out     = get_energy_flux(1,1)

     ! Pi_e/Pi_GB
     glf23_elec_mflux_out = get_stress_tor(1,1)

     ! S_e/S_GB
     glf23_elec_expwd_out = get_exchange(1,1)

     ! Ions

     do i_ion=1,2

        ! Gammai/Gamma_GB
        glf23_ion_pflux_out(i_ion) = get_particle_flux(i_ion+1,1)

        ! Qi/Q_GB
        glf23_ion_eflux_out(i_ion)     = get_energy_flux(i_ion+1,1)

        ! Pi_i/Pi_GB
        glf23_ion_mflux_out(i_ion) = get_stress_tor(i_ion+1,1)

        ! S_i/S_GB
        glf23_ion_expwd_out(i_ion) = get_exchange(i_ion+1,1)

     enddo

  else

     ! Run single-ky linear stability
     call glf2d

     ! Collect linear eigenvalues
     do n=1,2
        glf23_eigenvalue_out(n) = get_frequency(n) + xi*get_growthrate(n)
     enddo

  endif

!  call get_error_status(glf23_error_message,glf23_error_status)

end subroutine glf23_run
