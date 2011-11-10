
subroutine vgen_compute_neo(i,vtor_diff, rotation_model, er0, omega, omega_deriv)

  use vgen_globals
  use neo_interface
  use EXPRO_interface

  implicit none

  integer, intent(in) :: i
  integer, intent(in) :: rotation_model
  real, intent(in)    :: er0            ! kV
  real, intent(in)    :: omega          ! 1/s
  real, intent(in)    :: omega_deriv    ! 1/(m s)
  real, intent(out)   :: vtor_diff      ! vtor_exp - vtor_neo (m/s) 

  integer :: j, n
  real :: cc, loglam

  ! Set the local NEO input parameters
  
  neo_silent_flag_in    = 1

  ! Normalizations
  
  dens_norm = EXPRO_ne(i)
  temp_norm = EXPRO_te(i)
  mass_norm = neo_mass_1_in
  vth_norm  = sqrt(temp_norm * temp_norm_fac &
       / (mass_norm * mass_deuterium)) &
       * 1.0e4 / EXPRO_rmin(EXPRO_n_exp)

  cc = sqrt(2.0) * pi * charge_norm_fac**4 &
       * 1.0 / (4.0 * pi * 8.8542)**2 &
       * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
       * 1e9
  
  ! Rotation paramters
  
  neo_rotation_model_in  = rotation_model
  neo_dphi0dr_in         = -er0 / EXPRO_grad_r0(i) &
       * EXPRO_rmin(EXPRO_n_exp) / temp_norm
  neo_omega_rot_in       = omega / vth_norm 
  neo_omega_rot_deriv_in = omega_deriv / vth_norm * EXPRO_rmin(EXPRO_n_exp) 
  
  ! Geometry parameters
  
  neo_rmin_over_a_in = EXPRO_rmin(i)/EXPRO_rmin(EXPRO_n_exp)
  neo_rmaj_over_a_in = EXPRO_rmaj(i)/EXPRO_rmin(EXPRO_n_exp)
  neo_q_in           = abs(EXPRO_q(i))
  neo_shear_in       = EXPRO_s(i)
  neo_shift_in       = EXPRO_drmaj(i)
  neo_kappa_in       = EXPRO_kappa(i)
  neo_s_kappa_in     = EXPRO_skappa(i)
  neo_delta_in       = EXPRO_delta(i)
  neo_s_delta_in     = EXPRO_sdelta(i)
  neo_zeta_in        = EXPRO_zeta(i) 
  neo_s_zeta_in      = EXPRO_szeta(i) 
  neo_zmag_over_a_in = EXPRO_zmag(i)/EXPRO_rmin(EXPRO_n_exp)
  neo_s_zmag_in      = EXPRO_dzmag(i) 
  
  neo_rho_star_in = sqrt(temp_norm * temp_norm_fac &
       * mass_norm * mass_deuterium) &
       / (charge_norm_fac * abs(EXPRO_bunit(i))) &
       * 1.0e-4 / EXPRO_rmin(EXPRO_n_exp)
  
  if(neo_equilibrium_model_in == 3) then
     neo_geo_ny_in = EXPRO_nfourier
     neo_geo_yin_in(:,:) = 0.0
     do n=0,EXPRO_nfourier
        do j=1,4  
           neo_geo_yin_in(j,n)   = EXPRO_geo(j,n,i)/EXPRO_rmin(EXPRO_n_exp)
           neo_geo_yin_in(j+4,n) = EXPRO_dgeo(j,n,i)
        enddo
     enddo
  endif
  
  ! Species-dependent parameters
  
  neo_ne_ade_in = EXPRO_ne(i) / dens_norm
  neo_te_ade_in = EXPRO_te(i) / temp_norm
  
  loglam = 24.0 - log(sqrt(EXPRO_ne(i)*1e13)/(EXPRO_te(i)*1000))
  
  ! Species 1 -- must be an ion species
  neo_dens_1_in = EXPRO_ni_new(i) / dens_norm
  neo_dlnndr_1_in = EXPRO_dlnnidr_new(i) * EXPRO_rmin(EXPRO_n_exp)
  neo_temp_1_in = EXPRO_ti(1,i) / temp_norm
  neo_dlntdr_1_in = EXPRO_dlntidr(1,i) * EXPRO_rmin(EXPRO_n_exp)
  neo_nu_1_in = cc * loglam * neo_dens_1_in * dens_norm * neo_z_1_in**4 &
       / (sqrt(neo_mass_1_in) * (neo_temp_1_in * temp_norm)**1.5) &
       / vth_norm
  
  ! Species 2
  
  if(neo_n_species_in >= 2) then
     if(neo_z_2_in /= -1) then
        neo_dens_2_in = EXPRO_ni(2,i) / dens_norm
        neo_dlnndr_2_in = EXPRO_dlnnidr(2,i) * EXPRO_rmin(EXPRO_n_exp)
        neo_temp_2_in = EXPRO_ti(1,i) / temp_norm ! assume equal ion temps
        neo_dlntdr_2_in = EXPRO_dlntidr(1,i) * EXPRO_rmin(EXPRO_n_exp)
     else
        neo_dens_2_in = EXPRO_ne(i) / dens_norm
        neo_temp_2_in = EXPRO_te(i) / temp_norm
        neo_dlnndr_2_in = EXPRO_dlnnedr(i) * EXPRO_rmin(EXPRO_n_exp)
        neo_dlntdr_2_in = EXPRO_dlntedr(i) * EXPRO_rmin(EXPRO_n_exp)
     endif
  endif
  
  ! Species 3
  
  if(neo_n_species_in >= 3) then
     if(neo_z_3_in /= -1) then
        neo_dens_3_in = EXPRO_ni(3,i) / dens_norm
        neo_dlnndr_3_in = EXPRO_dlnnidr(3,i) * EXPRO_rmin(EXPRO_n_exp)
        neo_temp_3_in = EXPRO_ti(1,i) / temp_norm ! assume equal ion temps
        neo_dlntdr_3_in = EXPRO_dlntidr(1,i) * EXPRO_rmin(EXPRO_n_exp)
     else
        neo_dens_3_in = EXPRO_ne(i) / dens_norm
        neo_temp_3_in = EXPRO_te(i) / temp_norm
        neo_dlnndr_3_in = EXPRO_dlnnedr(i) * EXPRO_rmin(EXPRO_n_exp)
        neo_dlntdr_3_in = EXPRO_dlntedr(i) * EXPRO_rmin(EXPRO_n_exp)
     endif
  endif
  
  ! Species 4
  
  if(neo_n_species_in >= 4) then
     if(neo_z_4_in /= -1) then
        neo_dens_4_in = EXPRO_ni(4,i) / dens_norm
        neo_dlnndr_4_in = EXPRO_dlnnidr(4,i) * EXPRO_rmin(EXPRO_n_exp)
        neo_temp_4_in = EXPRO_ti(1,i) / temp_norm  ! assume equal ion temps
        neo_dlntdr_4_in = EXPRO_dlntidr(1,i) * EXPRO_rmin(EXPRO_n_exp)
     else
        neo_dens_4_in = EXPRO_ne(i) / dens_norm
        neo_temp_4_in = EXPRO_te(i) / temp_norm
        neo_dlnndr_4_in = EXPRO_dlnnedr(i) * EXPRO_rmin(EXPRO_n_exp)
        neo_dlntdr_4_in = EXPRO_dlntedr(i) * EXPRO_rmin(EXPRO_n_exp)
     endif
  endif

  ! Species 5
  
  if(neo_n_species_in >= 5) then
     if(neo_z_5_in /= -1) then
        neo_dens_5_in = EXPRO_ni(5,i) / dens_norm
        neo_dlnndr_5_in = EXPRO_dlnnidr(5,i) * EXPRO_rmin(EXPRO_n_exp)
        neo_temp_5_in = EXPRO_ti(1,i) / temp_norm  ! assume equal ion temps
        neo_dlntdr_5_in = EXPRO_dlntidr(1,i) * EXPRO_rmin(EXPRO_n_exp)
     else
        neo_dens_5_in = EXPRO_ne(i) / dens_norm
        neo_temp_5_in = EXPRO_te(i) / temp_norm
        neo_dlnndr_5_in = EXPRO_dlnnedr(i) * EXPRO_rmin(EXPRO_n_exp)
        neo_dlntdr_5_in = EXPRO_dlntedr(i) * EXPRO_rmin(EXPRO_n_exp)
     endif
  endif

   ! Species 6
  
  if(neo_n_species_in >= 6) then
     if(neo_z_6_in /= -1) then
        neo_dens_6_in = EXPRO_ni(6,i) / dens_norm
        neo_dlnndr_6_in = EXPRO_dlnnidr(6,i) * EXPRO_rmin(EXPRO_n_exp)
        neo_temp_6_in = EXPRO_ti(1,i) / temp_norm  ! assume equal ion temps
        neo_dlntdr_6_in = EXPRO_dlntidr(1,i) * EXPRO_rmin(EXPRO_n_exp)
     else
        neo_dens_6_in = EXPRO_ne(i) / dens_norm
        neo_temp_6_in = EXPRO_te(i) / temp_norm
        neo_dlnndr_6_in = EXPRO_dlnnedr(i) * EXPRO_rmin(EXPRO_n_exp)
        neo_dlntdr_6_in = EXPRO_dlntedr(i) * EXPRO_rmin(EXPRO_n_exp)
     endif
  endif

  ! Run NEO
  call neo_run()

  vtor_diff = EXPRO_vtor(erspecies_indx,i) &
       - neo_vtor_dke_out(erspecies_indx) * vth_norm * EXPRO_rmin(EXPRO_n_exp)

end subroutine vgen_compute_neo
