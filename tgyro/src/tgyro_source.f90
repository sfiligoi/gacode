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

  integer :: i
  real, external :: sigv
  real, external :: sigvth


  !-------------------------------------------------------
  ! Source terms (erg/cm^3/s):
  !
  do i=1,n_r

     ! Alpha power
     !  - sigv in cm^3/s

     s_alpha(i) = (ni(1,i)/2)**2*sigv(ti(1,i)/1e3)*e_alpha

     ! Bremsstrahlung radiation
     ! - From NRL formulary 
     ! - 1 W/cm^3 = 1e7 erg/cm^3/s

     s_brem(i) = 1e7*1.69e-32*ne(i)**2*sqrt(te(i))*z_eff(i)

     ! Electron-ion energy exchange
     ! - Positive as defined on RHS of ion equation
     ! - Multiply formulary expression by (3/2)ne:

     s_exch(i) = 1.5*nu_exch(i)*ne(i)*k*(te(i)-ti(1,i))

     ! Auxiliary power

     s_aux(i) = (1.0-r(i)/r_min)**loc_aux_exp*loc_aux_amp

  enddo
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Powers in units of erg/s

  ! Get integrated alpha-power
  call tgyro_volume_int(s_alpha,p_alpha)

  ! Get integrated Bremsstrahlung power
  call tgyro_volume_int(s_brem,p_brem)

  ! Get integrated exchange power
  call tgyro_volume_int(s_exch,p_exch)

  ! Get integrated auxiliary power
  call tgyro_volume_int(s_aux,p_aux)
  !-------------------------------------------------------

  !-------------------------------------------------------
  select case (loc_scenario)

  case(1)

     ! Experimental, static exchange
     ! 
     ! Input power is fixed and exchange is purely 
     ! diagnostic

     p_exch(:) = p_exch_in(:)

     p_i(:) = p_i_in(:)
     p_e(:) = p_e_in(:)

  case(2)

     ! Experimental, dynamic exchange
     ! 
     ! Input power is fixed but exchange is computed 

     p_i(:) = p_i_in(:)+(p_exch(:)-p_exch_in(:))
     p_e(:) = p_e_in(:)-(p_exch(:)-p_exch_in(:))

  case(3)

     ! Reactor, with self-consistent power, radiation 
     ! and exchange; input auxiliary power.

     p_i = (1.0-loc_alpha_elec)*p_alpha & 
          +p_exch &                       
          +p_i_aux_in                       

     p_e = loc_alpha_elec*p_alpha  &
          -p_exch &      
          +p_e_aux_in &                   
          -p_brem                       

  case(4)

     ! Reactor, with self-consistent power, radiation 
     ! and exchange; model auxiliary power.

     p_i = (1.0-loc_alpha_elec)*p_alpha & 
          +p_exch &                       
          +0.9*p_aux                       

     p_e = loc_alpha_elec*p_alpha  &
          -p_exch &      
          +0.1*p_aux &                   
          -p_brem                       

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

  end select
  !------------------------------------------------

  !------------------------------------------------
  ! Target angular momentum fluxes in erg/cm^2
  ! 
  mflux_target(1) = 0.0 

  mflux_target(2:n_r) = 0.0
  !------------------------------------------------

  !------------------------------------------------
  ! Target fluxes in GB units
  eflux_i_target(:) = eflux_i_target(:)/q_gb(:)
  eflux_e_target(:) = eflux_e_target(:)/q_gb(:)
  pflux_e_target(:) = pflux_e_target(:)/gamma_gb(:)
  !------------------------------------------------

end subroutine tgyro_source
