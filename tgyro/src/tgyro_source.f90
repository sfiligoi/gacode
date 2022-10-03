!------------------------------------------------------------
! tgyro_source.f90
!
! PURPOSE:
!  Calculation of primitive source elements, and volume 
!  integration into required form.
!------------------------------------------------------------

subroutine tgyro_source 

  use tgyro_globals

  implicit none

  integer :: i_ion

  !-------------------------------------------------------
  ! Source terms (erg/cm^3/s):
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! 1. Alpha power
  call rad_alpha(ne,ni,te,ti,sn_alpha,s_alpha_i,s_alpha_e,frac_ai,e_cross,n_alpha,t_alpha,n_r,loc_n_ion)
  frac_ae = 1-frac_ai
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! 2. Bremsstrahlung and line radiation (s_brem,s_line)
  call rad_ion_adas(te,ne,ni,zi_vec,ion_name,s_brem,s_line,loc_n_ion,n_r)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! 3. Synchrotron radiation (s_sync) with reflection co.
  call rad_sync(aspect_rat,r_min,b_ref,ne,te,s_sync,n_r)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! 4. Classical electron-ion energy exchange
  !  - Positive as defined on RHS of ion equation
  !  - Multiply formulary expression by (3/2)ne:
  s_exch(:) = 1.5*nu_exch(:)*ne(:)*k*(te(:)-ti(1,:))
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! 5. Anomalous electron-ion energy exchange
  !  - Positive as defined on RHS of ion equation
  !  - Skip exchange with fast ions
  s_expwd(:) = 0.0
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        s_expwd(:) = s_expwd(:)+expwd_i_tur(i_ion,:)*s_gb(:)
     endif
  enddo
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Powers in units of erg/s

  ! Integrated alpha-power
  call tgyro_volume_int(s_alpha_i,p_i_fus)
  call tgyro_volume_int(s_alpha_e,p_e_fus)
  call tgyro_volume_int(sn_alpha,f_he_fus)

  ! Integrated Bremsstrahlung power
  call tgyro_volume_int(s_brem,p_brem)

  ! Integrated Synchrotron power
  call tgyro_volume_int(s_sync,p_sync)

  ! Integrated line power
  call tgyro_volume_int(s_line,p_line)

  ! Integrated collisional exchange power
  call tgyro_volume_int(s_exch,p_exch)

  ! Integrated anomalous exchange power
  call tgyro_volume_int(s_expwd,p_expwd)
  !-------------------------------------------------------

  !-------------------------------------------------------
  select case (loc_scenario)

  case (1)

     ! Experimental, static exchange

     p_i(:) = p_i_in(:) ! Total ion input power (including exchange)
     p_e(:) = p_e_in(:) ! Total electon input power (including exchange)

  case (2)

     ! Experimental, dynamic exchange

     p_i(:) = p_i_in(:) &              ! Total ion input power 
          +(p_exch(:)-p_exch_in(:)) &  ! Consistent e-i exchange
          +p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

     p_e(:) = p_e_in(:) &              ! Total electron input power 
          -(p_exch(:)-p_exch_in(:)) &  ! Consistent e-i exchange
          -p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

  case (3)

     ! Reactor with consistent alpha power, exchange, radiation and Ohmic heating

     p_i(:) = &
          +p_i_fus(:) &                ! Fusion power to ions
          +p_i_aux_in(:) &             ! Auxiliary ion heating [fixed]
          +p_exch(:) &                 ! Collisional exchange
          +p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

     p_e(:) = &
          +p_e_fus(:) &                ! Fusion power to electrons
          +p_e_aux_in(:) &             ! Auxiliary electron heating [fixed]
          +p_e_ohmic_in(:) &           ! Ohmic heating [fixed]
          -p_exch(:)   &               ! Collisional exchange
          -p_brem(:) &                 ! Bremsstrahlung radiation
          -p_sync(:) &                 ! Synchrotron radiation
          -p_line(:) &                 ! Line radiation 
          -p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

  end select
  !-------------------------------------------------------

  !------------------------------------------------
  ! Target energy fluxes in erg/s/cm^2 
  eflux_i_target(1) = 0.0 
  eflux_e_target(1) = 0.0

  eflux_i_target(2:n_r) = p_i(2:n_r)/volp(2:n_r) 
  eflux_e_target(2:n_r) = p_e(2:n_r)/volp(2:n_r)
  !------------------------------------------------

  !------------------------------------------------
  ! Target particle fluxes in 1/s/cm^2
  select case (loc_pflux_method)

  case (1)

     pflux_e_target(:) = 0.0

  case (2)

     pflux_e_target(1) = 0.0
     pflux_e_target(2:n_r) = f_b_in(2:n_r)/volp(2:n_r)

  case (3)

     pflux_e_target(1) = 0.0
     pflux_e_target(2:n_r) = (f_b_in(2:n_r) + f_w_in(2:n_r)) &
          /volp(2:n_r) 

  end select
  !------------------------------------------------

  !------------------------------------------------
  ! Target angular momentum fluxes in erg/cm^2
  ! 
  mflux_target(1) = 0.0 
  mflux_target(2:n_r) = mf_in(2:n_r)/volp(2:n_r)
  !------------------------------------------------

  !------------------------------------------------
  ! Target He ash flux in 1/s/cm^2
  !
  pflux_he_target(1) = 0.0
  pflux_he_target(2:n_r) = f_he_fus(2:n_r)/volp(2:n_r)
  !------------------------------------------------

  !------------------------------------------------
  ! Target fluxes in GB units
  eflux_i_target(:)  = eflux_i_target(:)/q_gb(:)
  eflux_e_target(:)  = eflux_e_target(:)/q_gb(:)
  pflux_e_target(:)  = pflux_e_target(:)/gamma_gb(:)
  mflux_target(:)    = mflux_target(:)/pi_gb(:)
  pflux_he_target(:) = pflux_he_target(:)/gamma_gb(:)
  !------------------------------------------------

end subroutine tgyro_source

