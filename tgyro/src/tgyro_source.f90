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
  real, external :: dtrate_dv
  real :: n_d,n_t
  real :: s_alpha
  integer, parameter :: iprint=0

  !-------------------------------------------------------
  ! Source terms (erg/cm^3/s):
  !
  do i=1,n_r

     ! Alpha power
     !  - sigv in cm^3/s

     if (tgyro_dt_method == 1) then
        ! Assume D and T given by ion 1 and ion 2 
        ! (order doesn't matter)
        n_d = ni(1,i)
        n_t = ni(2,i)
     else
        ! Assume ion 1 is DT hybrid.
        n_d = 0.5*ni(1,i)
        n_t = 0.5*ni(1,i)
     endif

     s_alpha = n_d*n_t*sigv(ti(1,i)/1e3)*e_alpha
 !    s_alpha = n_d*n_t*dtrate_dv(ti(1,i)/1e3)*e_alpha

     s_alpha_i(i) = s_alpha*frac_ai(i)
     s_alpha_e(i) = s_alpha*frac_ae(i)

     ! Bremsstrahlung radiation
     ! - From NRL formulary 
     ! - 1 W/cm^3 = 1e7 erg/cm^3/s

     s_brem(i) = 1e7*1.69e-32*ne(i)**2*sqrt(te(i))*z_eff(i)

     ! Classical electron-ion energy exchange
     ! - Positive as defined on RHS of ion equation
     ! - Multiply formulary expression by (3/2)ne:

     s_exch(i) = 1.5*nu_exch(i)*ne(i)*k*(te(i)-ti(1,i))

     ! Anomalous electron-ion energy exchange
     ! - Positive as defined on RHS of ion equation
     ! - Average ion and -electron results to obtain best estimate  
     !   of ion value.

     s_expwd(i) = 0.5*(sum(expwd_i_tur(1:loc_n_ion,i))-expwd_e_tur(i))*s_gb(i)

  enddo
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Powers in units of erg/s

  ! Get integrated alpha-power
  call tgyro_volume_int(s_alpha_i,p_alpha_i)
  call tgyro_volume_int(s_alpha_e,p_alpha_e)

  if (i_proc_global == 0 .and. iprint == 1) then
     print '(17(1pe11.4,1x))',r(n_r)/r_min
     print '(17(1pe11.4,1x))',p_alpha_i(n_r)
     print '(17(1pe11.4,1x))',p_alpha_i_in(n_r)
     print *
     print '(17(1pe11.4,1x))',p_alpha_e(n_r)
     print '(17(1pe11.4,1x))',p_alpha_e_in(n_r)
     print *
     print '(17(1pe11.4,1x))',p_alpha_i(n_r)+p_alpha_e(n_r)
     print '(17(1pe11.4,1x))',p_alpha_i_in(n_r)+p_alpha_e_in(n_r)
  endif

  ! Get integrated Bremsstrahlung power
  call tgyro_volume_int(s_brem,p_brem)

  ! Get integrated collisional exchange power
  call tgyro_volume_int(s_exch,p_exch)

  ! Get integrated anomalous exchange power
  call tgyro_volume_int(s_expwd,p_expwd)
  !-------------------------------------------------------

  !-------------------------------------------------------
  select case (loc_scenario)

  case (1)

     ! Experimental, static exchange

     p_exch(:) = p_exch_in(:) ! Diagnostic exchange

     p_i(:) = p_i_in(:) ! Fixed total ion input power
     p_e(:) = p_e_in(:) ! Fixed total electon input power

  case (2)

     ! Experimental, dynamic exchange

     p_i(:) = p_i_in(:) &              ! Total ion input power 
          +(p_exch(:)-p_exch_in(:)) &  ! Consistent e-i exchange
          +p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

     p_e(:) = p_e_in(:) &              ! Total electron input power 
          -(p_exch(:)-p_exch_in(:)) &  ! Consistent e-i exchange
          -p_expwd(:)*tgyro_expwd_flag ! Turbulent exchange

  case (3)

     ! Reactor with consistent alpha power

     p_i(:) = p_i_in(:) &                   ! Total ion input power  
          +(p_exch(:)-p_exch_in(:)) &       ! Consistent e-i exchange
          +p_expwd(:)*tgyro_expwd_flag &    ! Turbulent exchange
          +(p_alpha_i(:)-p_alpha_i_in(:))   ! Consistent alpha power to ions

     p_e(:) = p_e_in(:) &                   ! Total electron input power
          -(p_exch(:)-p_exch_in(:)) &       ! Consistent e-i exchange
          -p_expwd(:)*tgyro_expwd_flag &    ! Turbulent exchange
          +(p_alpha_e(:)-p_alpha_e_in(:)) & ! Consistent alpha power to electrons
          -p_brem(:)                        ! Electron Bremsstrahlung

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
  ! Target fluxes in GB units
  eflux_i_target(:) = eflux_i_target(:)/q_gb(:)
  eflux_e_target(:) = eflux_e_target(:)/q_gb(:)
  pflux_e_target(:) = pflux_e_target(:)/gamma_gb(:)
  mflux_target(:)   = mflux_target(:)/pi_gb(:)
  !------------------------------------------------

end subroutine tgyro_source

real function dtrate_dv(ti)

  implicit none
  !\
  ! --------- 11/20/95 --- HSJ ---------------------------
  !  returns rate(cm**3/sec) of t(d,n)he4 reaction
  !  New Bosch & Hale rate coefficient:
  !  Bosch & Hale, Nuc. Fus., vol32, no.4 (1992) 611
  ! -------------------------------------------------------
  ! data for t(d,n)he4:
  !/
  real :: C1,C2,C3,C4,C5,C6,C7, B_gsq, mrcsq
  data  C1,C2,C3,C4,C5,C6,C7, B_gsq, mrcsq &
       / 1.17302e-9,1.51361e-2,7.51886e-2,4.60643e-3,1.3500e-2, &
       -1.06750e-4,1.36600e-5,1.182170e3,1124656.0/
  real :: theta,xsi
  real, intent(in) :: ti
  !\-------------------------------------------------
  ! neutrons produced by bulk plasma d-t fusion:
  !/
  theta  = ti*(C2+ti*(C4+ti*C6))
  theta  = theta/(1.0+ti*(C3+ti*(C5+ti*C7)))
  theta  = ti/(1.0-theta)
  xsi    = (B_gsq/(4.0*theta))**(1.0/3.0)
  dtrate_dv = C1*theta*sqrt(xsi/(mrcsq*ti**3))*exp(-3.0*xsi)
 
end function dtrate_dv
