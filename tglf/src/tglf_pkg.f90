      MODULE tglf_pkg
!
      IMPLICIT NONE
! output routines
      REAL,EXTERNAL :: get_growthrate,get_frequency
      REAL,EXTERNAL :: get_QL_particle_flux,get_QL_energy_flux
      REAL,EXTERNAL :: get_QL_stress_par,get_QL_stress_tor
      REAL,EXTERNAL :: get_QL_exchange
      REAL,EXTERNAL :: get_QL_phi,get_QL_density,get_QL_temperature
      REAL,EXTERNAL :: get_gaussian_width
      REAL,EXTERNAL :: get_R_unit
      REAL,EXTERNAL :: get_q_unit
      REAL,EXTERNAL :: get_wd_bar
      REAL,EXTERNAL :: get_b0_bar
      REAL,EXTERNAL :: get_ave_wd
      REAL,EXTERNAL :: get_ave_b0
      REAL,EXTERNAL :: get_particle_flux
      REAL,EXTERNAL :: get_energy_flux
      REAL,EXTERNAL :: get_stress_par
      REAL,EXTERNAL :: get_stress_tor
      REAL,EXTERNAL :: get_exchange
      REAL,EXTERNAL :: get_phi_bar
      REAL,EXTERNAL :: get_v_bar
      REAL,EXTERNAL :: get_n_bar
      REAL,EXTERNAL :: get_t_bar
      REAL,EXTERNAL :: get_n_bar_sum
      REAL,EXTERNAL :: get_t_bar_sum
      REAL,EXTERNAL :: get_Ne_Te_phase
      REAL,EXTERNAL :: get_phi_bar_sum
      REAL,EXTERNAL :: get_v_bar_sum
      REAL,EXTERNAL :: get_q_low
      REAL,EXTERNAL :: get_a_pol
      REAL,EXTERNAL :: get_a_tor
      REAL,EXTERNAL :: get_Bp0
      REAL,EXTERNAL :: get_R2_ave
      REAL,EXTERNAL :: get_B2_ave
      REAL,EXTERNAL :: get_RBt_ave
      REAL,EXTERNAL :: get_nky_out
      REAL,EXTERNAL :: get_eigenvalue_spectrum_out
      REAL,EXTERNAL :: get_field_spectrum_out
      REAL,EXTERNAL :: get_intensity_spectrum_out
      REAL,EXTERNAL :: get_flux_spectrum_out
      REAL,EXTERNAL :: get_ky_spectrum_out
      REAL,EXTERNAL :: get_DM
      REAL,EXTERNAL :: get_DR
      REAL,EXTERNAL :: get_DEP_parameters
      REAL,EXTERNAL :: get_wavefunction_out
!
      CONTAINS
!
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
      END MODULE tglf_pkg
