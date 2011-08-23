!-----------------------------------------------------------
! gyro_profile_init.f90
!
! PURPOSE:
!  Manage generation of radial profiles and associated 
!  quantities given radial_profile_method.
!
! NOTES:
!
!  radial_profile_method = 1: 
!    flat, profiles, s-alpha geometry
!  radial_profile_method = 5: 
!    flat profiles, s-alpha, model, or general geometry
!  radial_profile_method = 3: 
!    experimental profiles, model, or general geometry
!
!  r_e is the uniform REFERENCE grid.
!  r is the (possibly nonuniform) PHYSICAL grid.
!-----------------------------------------------------------

subroutine gyro_profile_init

  use gyro_globals
  use gyro_profile_exp
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer :: ic 
  real :: loglam
  real :: cc
  real :: p_total
  !---------------------------------------------------

  !---------------------------------------------------
  ! Generate toroidal mode numbers:
  !
  ! Set for each processor subgroup as well.
  !
  ! Super:
  do in=1,n_n
     n(in) = n0+(in-1)*d_n
  enddo
  !
  ! Sub:
  in_1 = 1
  !
  in = i_group_1*n_n_1 + in_1
  n_1(in_1) = n(in)
  !
  ! For nonlinear convolution:
  n_max = n_n-1 
  !---------------------------------------------------

  !---------------------------------------------------
  ! Charge and mass for each species:
  !
  ! First, ions:
  !
  do is=1,n_spec-1
     z(is)  = z_vec(is)
     mu(is) = mu_vec(is) 
  enddo
  !
  ! Finally, electrons.  Before any profile mapping, these 
  ! have species index n_spec.
  !
  z(n_spec)  = z_vec(0)
  mu(n_spec) = mu_vec(0)
  !---------------------------------------------------

  !---------------------------------------------------
  ! Set species-specific energy_max or default to
  ! one energy_max for all species.
  !
  if (variable_egrid_flag == 1) then
     energy_max(1:n_ion) = energy_max_vec(1:n_ion)
     energy_max(indx_e) = energy_max_vec(0)
  else
     energy_max(:) = energy_max_vec(1)
  endif
  !---------------------------------------------------

  !---------------------------------------------------
  ! Defaults and initializations.
  !
  if (radial_profile_method /= 3) then

     ! For local runs, set some unused parameters:

     a_meters   = 1.0
     bt_exp     = 1.0
     n_grid_exp = 0

  endif
  !
  ! Location of "norm point" (index of central radius)
  !
  ir_norm = n_x/2+1-n_x_offset
  !
  ! Set some geometry parameters not defined for 
  ! circular equilibrium.  These will also be reset 
  ! when using radial_profile_method=3  
  !
  b_unit_s(:)  = -btccw
  kappa_s(:)   = kappa0
  delta_s(:)   = delta0
  zeta_s(:)    = zeta0
  s_kappa_s(:) = s_kappa0
  s_delta_s(:) = s_delta0
  s_zeta_s(:)  = s_zeta0
  rmaj_s(:)    = r_maj
  drmaj_s(:)   = drmaj0
  zmag_s(:)    = zmag0
  dzmag_s(:)   = dzmag0

  !-------------------------------------------------------------------
  ! General geometry Fourier coefficients
  !
  if (num_equil_flag == 1) then
     do i=1,n_x
        a_fourier_geo_s(:,:,i) = a_fourier_geo(:,:)
     enddo
  else
     if (udsymmetry_flag == 1) then
        call send_message('INFO: Forcing up-down symmetry (UDSYMMETRY_FLAG=1).')
        zmag_s(:)  = 0.0
        dzmag_s(:) = 0.0
     endif
  endif
  !-------------------------------------------------------------------

  beta_unit_s(:) = betae_unit*(sum(n_vec(0:n_spec-1)*t_vec(0:n_spec-1)))

  csda_s(:)   = sqrt(1.6022e-16/(2.0*kg_proton))
  z_eff_s (:) = z_eff
  mach_s(:)   = 0.0
  w0_s(:)     = 0.0
  w0p_s(:)    = 0.0
  nu_s(:,:)   = 1.0
  !---------------------------------------------------

  !---------------------------------------------------
  ! Determine radial box size, etc.
  !
  call gyro_radial_simulation_box
  !---------------------------------------------------

  rhosda_s(:) = rho_star

  select case (radial_profile_method)

  case (1,5)

     !-----------------------------
     ! standard flat profiles (no variation)
     !-----------------------------

     ! Explicit sign convention
     q0 = q0*(ipccw)*(btccw)

     do i=1,n_x

        shat_s(i) = s0

        q(i) = q0+q0*s0/r0*(r(i)-r0)

        ! Get correct helicity

        r_s(i) = r0
        q_s(i) = q0

        rmaj_s(i) = r_maj

        ! Ions

        do is=1,n_ion
           den_s(is,i)    = n_vec(is)
           tem_s(is,i)    = t_vec(is)
           dlnndr_s(is,i) = dlnndr_vec(is) 
           dlntdr_s(is,i) = dlntdr_vec(is)
        enddo

        ! Electrons

        den_s(n_spec,i)    = n_vec(0)
        tem_s(n_spec,i)    = t_vec(0)
        dlnndr_s(n_spec,i) = dlnndr_vec(0)
        dlntdr_s(n_spec,i) = dlntdr_vec(0)

        ! Rotation with correct helicity

        gamma_p_s(i) = pgamma0
        gamma_e_s(i) = gamma_e
        mach_s(i)    = mach0

     enddo

  case (3) 

     !----------------------------------------
     ! experimental profiles and real geometry
     !---------------------------------------- 

     !---------------------------------------------------------------
     ! Read/broadcast experimental profiles on profile (_p) grid 
     ! and interpolate to simulation (_s) grid.  In general, the 
     ! _s grid is much finer than the _p grid.
     !
     ! At this point r = r_e and dr_eodr = 1.0.
     !
     call gyro_read_experimental_profiles
     !---------------------------------------------------------------

     call gyro_map_experimental_profiles

     !----------------------------------------------------------
     ! Recompute box length and gridpoints to get periodic case 
     ! right when using experimental profiles:
     !
     if (boundary_method == 1) then

        call send_message('INFO: Remapping into periodic domain.')

        if (box_multiplier < 0.0) then
           box_multiplier = 1.0
           ! rho_star < 0 if btccw = 1.
           rho_star = abs(l_y*r0/(n_ref*q_s(ir_norm)))*(-btccw)
           rhosda_s(:) = rho_star
        endif

        x_length = abs(box_multiplier*r0/(n_ref*q_s(ir_norm)*shat_s(ir_norm)))

        ! Box ends at periodic point
        !
        ! x----x----x----x----p
        !           |
        !          r=r0

        d_x = x_length/n_x

        do i=1,n_x
           r(i) = r0-0.5*x_length+d_x*(i-1.0)
        enddo

        r_e(:)     = r(:)
        dr_eodr(:) = 1.0

     endif
     !------------------------------------------------------

     !------------------------------------------------------
     ! Need to set r_s before exiting make_profile_*.
     ! So, set it for nonperiodic case and reset later if 
     ! flat_profile_flag=1:
     r_s = r
     !------------------------------------------------------

     !------------------------------------------------------
     ! Modification of input temperature and gradients
     !
     ! Experimental chi_gb_(:) will increase by a factor of
     !  
     !                  1/(1-eps_lt_vec(:))
     !
     ! Gradients in simulation will decrease by a factor of 
     !
     !                    1-eps_lt_vec(:)
     !
     ! Similar results hold for eps_ln_vec.
     !
     do is=1,n_spec-1
        dlntdr_s(is,:) = dlntdr_s(is,:)*(1.0-eps_dlntdr_vec(is))
        dlnndr_s(is,:) = dlnndr_s(is,:)*(1.0-eps_dlnndr_vec(is))
        if (eps_dlntdr_vec(is) /= 0.0) then
           call send_message_real('INFO: Ti gradient RESCALED by: ',&
                1.0-eps_dlntdr_vec(is))
           call logint(tem_s(is,:),dlntdr_s(is,:),r_s,n_x,ir_norm)
           call send_message('INFO: Ti profile REINTEGRATRED')
        endif
        if (eps_dlnndr_vec(is) /= 0.0) then
           call send_message_real('INFO: ni gradient RESCALED by: ',&
                1.0-eps_dlnndr_vec(is))
           call logint(den_s(is,:),dlntdr_s(is,:),r_s,n_x,ir_norm)
           call send_message('INFO: ni profile REINTEGRATRED')
        endif
     enddo
     !
     dlntdr_s(n_spec,:) = dlntdr_s(n_spec,:)*(1.0-eps_dlntdr_vec(0))
     dlnndr_s(n_spec,:) = dlnndr_s(n_spec,:)*(1.0-eps_dlnndr_vec(0))
     if (eps_dlntdr_vec(0) /= 0.0) then
        call send_message_real('INFO: Te gradient RESCALED by: ',&
             1.0-eps_dlntdr_vec(0))

        call logint(tem_s(n_spec,:),dlntdr_s(n_spec,:),r_s,n_x,ir_norm)
        call send_message('INFO: Te profile REINTEGRATRED')
     endif
     if (eps_dlnndr_vec(0) /= 0.0) then
        call send_message_real('INFO: ne gradient RESCALED by: ',&
             1.0-eps_dlnndr_vec(0))
        call logint(den_s(n_spec,:),dlnndr_s(n_spec,:),r_s,n_x,ir_norm)
        call send_message('INFO: ne profile REINTEGRATRED')
     endif
     !
     if ((sum(abs(eps_dlntdr_vec(:))) + sum(abs(eps_dlnndr_vec(:)))) /= 0.0) then
        call send_message('INFO: profiles changed, recalculating beta_unit')
        do i=1,n_x
           p_total = sum(den_s(:,i)*tem_s(:,i))
           beta_unit_s(i) = 400.0*p_total/(1e5*b_unit_s(i)**2)
        enddo
     endif
     !------------------------------------------------------

  end select ! profile_method

  !--------------------------------------------------------------------
  ! Overwrite Z_eff with self-consistent value if z_eff_method set:
  !
  if (z_eff_method == 2) then

     z_eff_s(:) = 0.0

     do is=1,n_spec-1
        z_eff_s(:) = z_eff_s(:)+z(is)**2*den_s(is,:)/den_s(n_spec,:)
     enddo

  endif
  !--------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! When electron_method=4, n_spec is naturally overlayed on n_ion, but 
  ! for electron_method=3, we must permute species
  !
  if (electron_method == 3) then
     z(1:n_spec)          = z(n_spec:1:-1)
     mu(1:n_spec)         = mu(n_spec:1:-1)
     den_s(1:n_spec,:)    = den_s(n_spec:1:-1,:)
     tem_s(1:n_spec,:)    = tem_s(n_spec:1:-1,:)
     dlntdr_s(1:n_spec,:) = dlntdr_s(n_spec:1:-1,:)
     dlnndr_s(1:n_spec,:) = dlnndr_s(n_spec:1:-1,:)
  endif
  !------------------------------------------------------------------

  !-------------------------------------------------------------
  ! Values of functions at the norm point.
  !
  ! These are reset from the intialized values in 
  ! make_profiles.
  !
  r_norm     = r(ir_norm)
  q_norm     = q_s(ir_norm)
  shat_norm  = shat_s(ir_norm)
  rhos_norm  = rhosda_s(ir_norm)
  csda_norm  = csda_s(ir_norm)

  tem_norm = tem_s(indx_e,ir_norm)
  den_norm = den_s(indx_e,ir_norm)

  b_unit_norm = b_unit_s(ir_norm)
  r_maj_norm  = rmaj_s(ir_norm)
  !
  ! betae_unit_norm enters Ampere's law
  !
  pr_s(:,:) = den_s(:,:)*tem_s(:,:)
  !
  betae_unit_norm = beta_unit_s(ir_norm)*pr_s(indx_e,ir_norm) &
       /sum(pr_s(1:n_spec,ir_norm))
  !-------------------------------------------------------------

  !------------------------------------------------------
  ! FINAL NORMALIZATIONS:
  !
  ! Convert b_unit_s, tem_s and den_s to units of   
  ! central b_unit, Te and ne.
  !
  b_unit_s(:) = b_unit_s(:)/b_unit_norm
  tem_s(:,:)  = tem_s(:,:)/tem_norm
  den_s(:,:)  = den_s(:,:)/den_norm
  !
  alpha_s(:,:) = den_s(:,:)/tem_s(:,:)
  !------------------------------------------------------

  !------------------------------------------------------
  ! Total pressure gradient
  !
  do i=1,n_x
     dlnpdr_s(i) = sum(pr_s(1:n_spec,i)*(dlnndr_s(1:n_spec,i)+&
          dlntdr_s(1:n_spec,i)))/sum(pr_s(1:n_spec,i))
  enddo
  !
  ! beta_star = -(8 pi)/(B_unit**2) dp/dr
  !
  if (geo_fastionbeta_flag == 0) then

     ! Pressure from species sum

     beta_star_s(:) = beta_unit_s(:)*dlnpdr_s(:)*geo_betaprime_scale         
     if (geo_betaprime_scale /= 1.0) then
        call send_message_real(&
             'INFO: Scaling dp/dr in GEO by: ',geo_betaprime_scale)
     endif

  else

     ! Pressure from total pressure including fast ions

     beta_star_s(:) = beta_unit_ptot_s(:)*dlnptotdr_s(:)
     call send_message('INFO: Using total dp/dr (+ fast ions) in GEO.')

  endif
  !------------------------------------------------------

  !----------------------------------------------------------
  ! Sign manipulations: up to now, we k*rho can have either 
  ! sign depending on field and current orientation.
  ! Now, we manipulate the sign of n to make k*rho >=0:
  !
  if (q_s(ir_norm)*rhos_norm < 0) then
     n_1(:) = -n_1(:)
     n0     = -n0
     n(:)   = -n(:)
  endif
  !
  ! Define krho_i = k_theta * rho_{s,unit} 
  !
  do i=1,n_x
     krho_i(:,i) = n_1(:)*q_s(i)/r_s(i)*rhos_norm/b_unit_s(i)
  enddo
  !----------------------------------------------------------

  !---------------------------------------------------------
  ! Rotation parameter scaling and definition of omega_eb_s
  ! 
  !  omega0 = - c (dphi/dpsi) 
  !
  !       M = omega_0 R_0/c_s
  ! gamma_p = -R_0 d(omega_0)/dr
  ! gamma_e = -(r/q) d(omega_0)/dr
  !
  if (radial_profile_method == 3) then

     mach_s(:)     = mach_s(:)*mach0_scale
     gamma_e_s(:)  = gamma_e_s(:)/csda_norm*doppler_scale
     gamma_p_s(:)  = gamma_p_s(:)/csda_norm*pgamma0_scale
     w0_s(:)       = w0_s(:)/csda_norm*doppler_scale

     omega_eb_s(:) = -n_1(in_1)*w0_s(:)

     if (mach0_scale /= 1.0) then
        call send_message_real('INFO: Scaling experimental Mach number by ',&
             mach0_scale)
     endif

     if (pgamma0_scale /= 1.0) then
        call send_message_real('INFO: Scaling experimental gamma_p by ',&
             pgamma0_scale)
     endif

     if (doppler_scale /= 1.0) then
        call send_message_real('INFO: Scaling experimental gamma_e by ',&
             doppler_scale)
     endif

     !---------------------------------------------------------
     ! Subtraction of a constant rotational velocity 
     omega_eb_s(:) = omega_eb_s(:)-omega_eb_s(ir_norm)
     !---------------------------------------------------------

  else

     omega_eb_s(:) = gamma_e_s(ir_norm)*n_1(in_1)*q_s(:)/r_norm*(r(:)-r_norm)

  endif
  !---------------------------------------------------------

  !----------------------------------------------------------
  ! Complete definition of profile functions:
  !
  if (flat_profile_flag == 0 .and. radial_profile_method == 3) then

     ! Profiles have RADIAL VARIATION

     do i=1,n_x

        ! Normal definitions

        q(i) = q_s(i)

     enddo

  else

     ! Profiles are FLAT

     do i=1,n_x

        if (radial_profile_method == 3) then

           r_s(i)     = r_norm
           rmaj_s(i) = rmaj_s(ir_norm)

           ! q, as always, needs a linear form here:

           q(i) = q_s(ir_norm)+   &
                q_s(ir_norm)/r(ir_norm)*shat_s(ir_norm)*(r(i)-r(ir_norm))

        endif

        rhogrid_s(i)    = rhogrid_s(ir_norm)
        q_s(i)          = q_s(ir_norm)
        kappa_s(i)      = kappa_s(ir_norm)
        delta_s(i)      = delta_s(ir_norm)
        zeta_s(i)       = zeta_s(ir_norm)
        z_eff_s(i)      = z_eff_s(ir_norm)
        b_unit_s(i)     = b_unit_s(ir_norm)
        rhosda_s(i)     = rhosda_s(ir_norm)
        csda_s(i)       = csda_s(ir_norm)
        drmaj_s(i)      = drmaj_s(ir_norm)
        s_delta_s(i)    = s_delta_s(ir_norm)
        s_zeta_s(i)     = s_zeta_s(ir_norm)
        s_kappa_s(i)    = s_kappa_s(ir_norm)
        shat_s(i)       = shat_s(ir_norm)
        beta_unit_s(i)  = beta_unit_s(ir_norm)
        beta_star_s(i)  = beta_star_s(ir_norm)

        gamma_e_s(i) = gamma_e_s(ir_norm)
        gamma_p_s(i) =  gamma_p_s(ir_norm)
        mach_s(i)    = mach_s(ir_norm)

        tem_s(:,i)    = tem_s(:,ir_norm)
        den_s(:,i)    = den_s(:,ir_norm)
        dlntdr_s(:,i) = dlntdr_s(:,ir_norm)
        dlnndr_s(:,i) = dlnndr_s(:,ir_norm)
        dlnpdr_s(i)   = dlnpdr_s(ir_norm)
        pr_s(:,i)     = pr_s(:,ir_norm)
        alpha_s(:,i)  = alpha_s(:,ir_norm)
        krho_i(:,i)   = krho_i(:,ir_norm)

        omega_eb_s(:) = gamma_e_s(ir_norm)*&
             n_1(in_1)*q_norm/r_norm*(r(:)-r_norm)

        a_fourier_geo_s(:,:,i) = a_fourier_geo_s(:,:,ir_norm)

     enddo

  endif
  !----------------------------------------------------------

  !------------------------------------------------------
  ! Compute true Debye length at box center if profile
  ! data is available:
  !
  if (radial_profile_method == 3) then

     !----------------------------------------------------------
     ! Debye length (from NRL plasma formulary):
     !
     lambda_debye = lambda_debye_scale*7.43* &
          sqrt((1e3*tem_norm)/(1e13*den_norm))/a_meters
     !----------------------------------------------------------

  endif
  !------------------------------------------------------

  !-------------------------------------------------------------
  ! Compute base collision rates
  !
  !          4 pi ne e^4 ne ln(Lambda)
  !  nu_ei = ------------------------- 
  !            me^(1/2) (2 Te)^(3/2)
  !
  !            nu_ei
  !  nu_e = ----------- { Zeff + H[sqrt(enhat_e)] } -> nu_s(indx_e,:)
  !         enhat^(3/2)

  if (radial_profile_method == 3) then

     ! Numerical coefficient
     cc = 4*pi*(4.8032e-10)**4/(2*1.6022e-12)**1.5/(9.1094e-28)**0.5  

     ! Coulomb logarithm
     loglam = 24.0-log(sqrt(den_norm*1e13)/(tem_norm*1e3))

     do is=1,n_spec

        ! Collision rate (1/sec)
        nu_s(is,:) = cc*loglam*(den_norm*den_s(is,:)*1e13)/&
             (tem_norm*tem_s(is,:)*1e3)**1.5*&
             mu(is)/mu(indx_e)*z(is)**2

        ! Express in dimensionless GYRO units:
        nu_s(is,:) = nu_s(is,:)/csda_norm

     enddo

  else

     do is=1,n_spec

        ! Collision rate (GYRO units):  
        nu_s(is,:) = nu_ei*(den_s(is,:)/den_s(indx_e,:))/ &
             (tem_s(is,:)/tem_s(indx_e,:))**1.5*&
             mu(is)/mu(indx_e)*z(is)**2

     enddo

  endif
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Apply collisional scaling factors, and determine number of 
  ! collisional species:
  !
  if (nu_ii_scale /= 1.0) then
     call send_message_real('INFO: nu_ii rescaled by: ',nu_ii_scale)
  endif
  if (nu_ei_scale /= 1.0) then
     call send_message_real('INFO: nu_ei rescaled by: ',nu_ei_scale)
  endif
  !
  ic = 0
  c_map(:) = 0.0
  do is=1,n_spec

     if (is /= indx_e) then 
        nu_s(is,:) = nu_s(is,:)*nu_ii_scale
     else
        nu_s(is,:) = nu_s(is,:)*nu_ei_scale
     endif

     if (nu_s(is,ir_norm) > 0.0) then
        ic = ic+1
        c_map(ic) = is
     endif

  enddo
  n_coll = ic
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Calculation of phase
  !
  do i=1,n_x
     phase(in_1,i) = exp(-pi_2*i_c*n_1(in_1)*q(i))
     angp(i) = n_1(in_1)*q(i)-int(n_1(in_1)*q(i))
     if (angp(i) > 0.5) angp(i) = angp(i)-1.0
  enddo
  !-------------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,"[gyro_profile_init done]"
  endif

end subroutine gyro_profile_init
