      MODULE tglf_pkg
!
      IMPLICIT NONE
      PRIVATE
! linear stability and transport model drivers
      PUBLIC :: tglf,tglf_TM
! gemetry setup routine
      PUBLIC :: tglf_setup_geometry
! input routines
      PUBLIC :: put_species,put_kys,put_signs
      PUBLIC :: put_gaussian_width,put_averages
      PUBLIC :: put_gradients,put_profile_shear
      PUBLIC :: put_switches,put_model_parameters
      PUBLIC :: put_s_alpha_geometry,put_Miller_geometry
      PUBLIC :: put_Fourier_geometry,put_ELITE_geometry
      PUBLIC :: put_eikonal, put_rare_switches
! output routines
      PUBLIC :: get_growthrate,get_frequency
      PUBLIC :: get_QL_particle_flux,get_QL_energy_flux
      PUBLIC :: get_QL_stress_par,get_QL_stress_tor 
      PUBLIC :: get_QL_exchange 
      PUBLIC :: get_QL_phi,get_QL_density,get_QL_temperature
      PUBLIC :: get_gaussian_width  
      PUBLIC :: get_R_unit
      PUBLIC :: get_q_unit
      PUBLIC :: get_wd_bar
      PUBLIC :: get_b0_bar
      PUBLIC :: get_ave_wd
      PUBLIC :: get_ave_b0
      PUBLIC :: get_particle_flux
      PUBLIC :: get_energy_flux
      PUBLIC :: get_stress_par
      PUBLIC :: get_stress_tor
      PUBLIC :: get_exchange
      PUBLIC :: get_phi_bar
      PUBLIC :: get_v_bar
      PUBLIC :: get_n_bar
      PUBLIC :: get_t_bar
      PUBLIC :: get_n_bar_sum
      PUBLIC :: get_t_bar_sum
      PUBLIC :: get_Ne_Te_phase
      PUBLIC :: get_phi_bar_sum
      PUBLIC :: get_v_bar_sum
      PUBLIC :: get_q_low
      PUBLIC :: get_a_pol
      PUBLIC :: get_a_tor
      PUBLIC :: get_Bp0
      PUBLIC :: get_R2_ave
      PUBLIC :: get_B2_ave
      PUBLIC :: get_RBt_ave
      PUBLIC :: get_DM
      PUBLIC :: get_DR
      PUBLIC :: get_wavefunction_out
      PUBLIC :: write_tglf_input
      PUBLIC :: write_wavefunction_out
!
      CONTAINS
!
      SUBROUTINE tglf
      USE tglf_global
      IMPLICIT NONE
      INTEGER :: i
!
      do i=1,7
        trace_path(i)=0
      enddo
!
      if(find_width_in)then
        trace_path(1)=1
        call tglf_max
      else
        call tglf_LS
      endif
      find_width_in=.FALSE.
!
      END SUBROUTINE tglf
!
      include 'tglf_inout.f90'
!      include 'tglf_geometry.f90'
!      include 'tglf_matrix.f90'
!      include 'tglf_LS.f90'
!      include 'tglf_max.f90'
      include 'tglf_setup_geometry.f90'
      include 'tglf_TM.f90'
!
      END MODULE tglf_pkg
