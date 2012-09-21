subroutine tgyro_profile_functions

  use tgyro_globals

  implicit none

  integer :: i_ion
  real :: c_exch
  real, dimension(n_r) :: loglam


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
     w0p(:) = f_rot(:)*w0p_norm
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
  gamma_gb(:) = ne(:)*c_s(:)*(rho_s(:)/r_min)**2

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

  ! NOTE: 
  ! c_exch = 1.8e-19 is the formulary exch. coefficient
  c_exch = 2.0*(4.0/3)*sqrt(2.0*pi)*e**4/k**1.5

  ! nu_exch in 1/s
  nu_exch(:) = 0.0
  do i_ion=1,loc_n_ion
     nu_exch(:) = nu_exch(:)+c_exch*sqrt(me*mi(i_ion))*zi_vec(i_ion)**2 &
          *ni(i_ion,:)*loglam(:)/(me*ti(i_ion,:)+mi(i_ion)*te(:))**1.5
  enddo

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

end subroutine tgyro_profile_functions
