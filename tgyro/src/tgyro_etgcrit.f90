subroutine tgyro_etgcrit

  use tgyro_globals

  implicit none

  ! Input parameters
  chi_gb_in = chi_gb(i_r)
  etae_in   = dlntedr(i_r)/dlnnedr(i_r)
  alte_in   = r_min*dlnnedr(i_r)

  ! Diffusivity model
  chi_e_out = 1.5*(etae_in-1.4)*(alte_in/60.0)*chi_gb_in
  d_e_out   = 0.02*chi_e_out

  ! Conversion to output fluxes
  etgcrit_elec_eflux_out = ne(i_r)*te(i_r)*dlntedr(i_r)*chi_e_out
  etgcrit_elec_pflux_out = ne(i_r)*dlnnedr(i_r)*d_e_out
  
  etgcrit_ion_eflux_out(1:loc_n_ion) = 0.0
  etgcrit_ion_pflux_out(1:loc_n_ion) = 0.0
  
end subroutine tgyro_etgcrit
