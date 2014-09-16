subroutine tgyro_profile_functions

  use tgyro_globals

  implicit none

  integer :: i
  integer :: i_ion
  real :: c_exch
  real, dimension(n_r) :: loglam
  real, dimension(n_r) :: c_a
  real, dimension(n_r) :: x_a
  real, external :: sivukhin

  ! Note flag to only evolve only gradients
  if (loc_evolve_grad_only_flag == 0 .and. &
       (loc_lock_profile_flag == 0 .or. i_tran > 0)) then

     !-------------------------------------------
     ! Integrate gradients to obtain profiles:
     !
     do i_ion=1,loc_n_ion
        ! ni in 1/cm^3
        call logint(ni(i_ion,:),dlnnidr(i_ion,:),r,n_r,i_bc)
        ! ti in eV
        call logint(ti(i_ion,:),dlntidr(i_ion,:),r,n_r,i_bc)
     enddo
     !
     ! ne in 1/cm^3
     call logint(ne,dlnnedr,r,n_r,i_bc)
     !
     ! te in eV
     call logint(te,dlntedr,r,n_r,i_bc)

     ! w0 in rad/s
     if (loc_er_feedback_flag == 1) then
        select case(tgyro_er_bc)
        case (1)
           ! Fixed derivative of w0 
           f_rot(1) = f_rot(1)
        case (2)
           ! Zero derivative of w0, similar to profiles
           f_rot(1) = 0.0
        case(3)
           ! Zero second derivative
           f_rot(1) = f_rot(2)
        end select
     endif

     ! NOTE (see tgyro_init_profiles)
     ! f_rot [1/cm] = w0p/w0_norm
     !
     ! w0_norm = c_s/R_maj at r=0.

     w0p(:) = f_rot(:)*w0_norm
     call linint(w0,w0p,r,n_r,i_bc)
     !-------------------------------------------

  endif

  ! Thermal velocities in cm/s
  v_i(:) = sqrt(k*ti(1,:)/mi(1))
  c_s(:) = sqrt(k*te(:)/mi(1))

  ! Thermal gyroradii in cm
  rho_i(:) = v_i(:)/(e*b_unit(:)/(mi(1)*c)) 
  rho_s(:) = c_s(:)/(e*b_unit(:)/(mi(1)*c))

  ! Gyrobohm unit diffusivity (cm^2/s)
  chi_gb(:) = rho_s(:)**2*c_s(:)/r_min

  ! Gyrobohm unit particle flux (1/cm^2/s)
  gamma_gb(:)  = ne(:)*c_s(:)*(rho_s(:)/r_min)**2

  ! Gyrobohm unit momentum flux (erg/cm^2)
  pi_gb(:) = ne(:)*k*te(:)*r_min*(rho_s(:)/r_min)**2

  ! Gyrobohm unit energy flux (erg/cm^2/s)
  q_gb(:) = ne(:)*k*te(:)*c_s(:)*(rho_s(:)/r_min)**2

  ! Gyrobohm unit exchange power density (erg/cm^3/s)
  s_gb(:) = ne(:)*k*te(:)*(c_s(:)/r_min)*(rho_s(:)/r_min)**2

  ! Coulomb logarithm
  loglam(:) = 24.0-log(sqrt(ne(:))/te(:))

  ! 1/tau_ii (Belli 2008) in 1/s
  do i_ion=1,loc_n_ion
     nui(i_ion,:) = sqrt(2.0)*pi*ni(i_ion,:)*(zi_vec(i_ion)*e)**4*loglam(:) &
          /(sqrt(mi(i_ion))*(k*ti(i_ion,:))**1.5)
  enddo

  ! 1/tau_ee (Belli 2008) in 1/s
  nue(:) = sqrt(2.0)*pi*ne(:)*e**4*loglam(:) &
       /(sqrt(me)*(k*te(:))**1.5)

  ! Hinton-Hazeltine scattering rates (one ion):
  nui_HH(:) = 4.0/(3*sqrt(2.0*pi))*nui(1,:)
  nue_HH(:) = 4.0/(3*sqrt(pi))*nue(:)*(ni(1,:)*zi_vec(1)**2/ne(:))

  ! INVERSE of nue_star
  nue_star(:) = (r(:)/r_maj(:))**1.5/abs(q(:))/nue_HH(:)* &
       (c_s(:)*sqrt(mi(1)/me))/r_maj(:)/z_eff(:)

  ! NOTE: 
  ! c_exch = 1.8e-19 is the formulary exch. coefficient
  c_exch = 2.0*(4.0/3)*sqrt(2.0*pi)*e**4/k**1.5

  ! nu_exch in 1/s
  nu_exch(:) = 0.0
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        nu_exch(:) = nu_exch(:)+c_exch*sqrt(me*mi(i_ion))*zi_vec(i_ion)**2 &
             *ni(i_ion,:)*loglam(:)/(me*ti(i_ion,:)+mi(i_ion)*te(:))**1.5
     endif
  enddo

  ! Alpha heating coefficients [Stix, Plasma Phys. 14 (1972) 367] 
  ! See in particular Eqs. 15 and 17.
  c_a(:) = 0.0
  do i_ion=1,loc_n_ion
     if (therm_flag(i_ion) == 1) then
        c_a(:) = c_a(:)+(ni(i_ion,:)/ne(:))*zi_vec(i_ion)**2/(mi(i_ion)/malpha)
     endif
  enddo
  e_cross(:) = k*te(:)*(4*sqrt(me/malpha)/(3*sqrt(pi)*c_a(:)))**(-2.0/3.0)
  x_a(:) = e_alpha/e_cross(:)
  do i=1,n_r
     frac_ai(i) = sivukhin(x_a(i))
  enddo
  frac_ae(:) = 1.0-frac_ai(:)

  ! Total pressure and beta (dimensionless) 
  pr(:) = ne(:)*k*te(:)
  do i_ion=1,loc_n_ion
     pr(:) = pr(:)+ni(i_ion,:)*k*ti(i_ion,:)
  enddo

  beta_unit(:)  = 8*pi*pr(:)/b_unit**2
  betae_unit(:) = beta_unit(:)*ne(:)*k*te(:)/pr(:)

  ! Pressure gradient inverse scale length (1/cm)
  dlnpdr(:) = ne(:)*k*te(:)*(dlnnedr(:)+dlntedr(:))/pr(:)
  do i_ion=1,loc_n_ion
     dlnpdr(:) = dlnpdr(:)+&
          ni(i_ion,:)*k*ti(i_ion,:)*(dlnnidr(i_ion,:)+dlntidr(i_ion,:))/pr(:)
  enddo

  !----------------------------------
  ! Functions connected with rotation
  !---------------------------------- 

  ! u00 (note that mach = u00/cs)
  u00(:) = r_maj(:)*w0(:)
  !
  ! gamma_eb (1/s)
  gamma_eb(:) = -r(:)/q(:)*w0p(:)
  !
  ! gamma_p (1/s)
  gamma_p(:)  = -r_maj(:)*w0p(:)

  !----------------------------------------------------------------------
  ! Trinity gyroBohm factor
  !
  q_tgb(:) = ni(1,:)*k*ti(1,:)*(sqrt(2.0)*v_i(:))*&
       (sqrt(2.0)*rho_i(:)*b_unit(:)/b_ref/r_min)**2
  !----------------------------------------------------------------------

end subroutine tgyro_profile_functions
