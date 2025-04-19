subroutine cgyro_check

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is
  character(len=1), dimension(7) :: ctag
  character(1000) :: outstr
  character(len=7) :: floatstr

  !-----------------------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_xi,2) /= 0) then 
     call cgyro_error('n_xi must be even.')
     return
  endif

  if (zf_test_mode == 0 .and. modulo(n_radial,box_size) /= 0) then 
     call cgyro_info('RESOLUTION WARNING -- n_radial not a multiple of box_size.')
  endif

  !if (zf_test_mode == 0 .and. n_radial < (n_toroidal-1)*box_size) then
  !   call cgyro_info('RESOLUTION WARNING -- n_radial < n*box_size.')
  !endif

  if (zf_test_mode == 0 .and. n_radial < box_size) then
     call cgyro_info('SEVERE RESOLUTION WARNING -- n_radial < box_size.')
     return
  endif

  if (modulo(n_theta,theta_plot) /= 0) then
     call cgyro_error('n_theta must be divisible by theta_plot.')
     return
  endif

  if (n_species > 11) then
     call cgyro_error('n_species cannot be greater than 11.')
     return
  endif

  if (h_print_flag == 1) then
     if (box_size > 1 .and. zf_test_mode == 0) then
        call cgyro_error('Distribution output not available for box_size > 1')
        return
     endif
  endif

  if (abs(px0) > 0.0) then
     if (box_size > 1) then
        call cgyro_error('Nonzero PXO not available for box_size > 1')
     else
        write(floatstr,'(f5.3)') 2*px0
        call cgyro_info('Finite-theta0 mode: theta0 = '//trim(floatstr)//' pi')
     endif
  endif
  !------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Time integration
  !
  select case (delta_t_method)
  case(0)
     call cgyro_info('Time integrator: RK4 4:4(3) [non-adaptive]')
  case(1)
     call cgyro_info('Time integrator: Cash-Karp 6:5(4) [adaptive]')
  case(2)
     call cgyro_info('Time integrator: Bogacki-Shampine 7:5(4) [adaptive]')
  case(3)
     call cgyro_info('Time integrator: Verner 10:7(6) [adaptive]')
  case default
     call cgyro_error('Invalid value for delta_t_method')
     return
  end select
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Electrons
  !
  select case (ae_flag)
  case(0)
     call cgyro_info('Using gyrokinetic electrons')
  case(1)
     call cgyro_info('Using adiabatic electrons')
  case default
     call cgyro_error('Invalid value for ae_flag')
     return
  end select
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Profile checks
  !
  select case (profile_model)

  case (1)
     call cgyro_info('Profile model: local input (input.cgyro)')

  case (2)
     call cgyro_info('Profile model: experimental (input.gacode)')

  case default
     call cgyro_error('Invalid value for profile_model')
     return

  end select
  !-----------------------------------------------------------------------

  if (profile_model == 2) then
     select case(quasineutral_flag)
     case(0)
        call cgyro_info('QN flag: Not enforcing quasi-neutrality')
     case(1)
        call cgyro_info('QN flag: Enforcing quasi-neutrality')
     case default
        call cgyro_error('Invalid value for quasineutral_flag')
        return
     end select
  endif

  !------------------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  

  case (1) 
     call cgyro_info('Equilibrium: s-alpha')
     if (profile_model == 2) then
        call cgyro_error('s-alpha equilibrium model not valid with experimental profiles')
        return
     endif

  case (2) 
     call cgyro_info('Equilibrium: Miller Extended Harmonic (MXH)')

  case (3) 
     call cgyro_info('Equilibrium: Fourier')

     if (geo_ny <= 0) then
        call cgyro_error('Fourier geometry coefficients missing.')
        return
     endif
     if (udsymmetry_flag == 1) then
        call cgyro_error('Cannot have UDSYMMETRY_FLAG=1 with general geometry.')
        return
     endif

  case default

     call cgyro_error('Invalid value for equilibrium_model')
     return

  end select
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Field consistency checks
  if (n_field > 1) then
     if (abs(betae_unit) < epsilon(0.0)) then
        call cgyro_error('BETAE_UNIT must be nonzero for electromagnetic simulation.')
        return
     endif
     if (ae_flag == 1) then
        call cgyro_error('Electrons must be gyrokinetic for electromagnetic simulation.')
        return
     endif
  endif
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Fields
  select case(n_field)

  case (1)
     call cgyro_info('Electrostatic fluctuations (Phi)')

  case (2)
     call cgyro_info('Transverse EM fluctuations (Phi,A_par)')

  case (3)
     call cgyro_info('Transverse and compressional EM fluctuations (Phi,A_par,B_par)')

  case default
     call cgyro_error('Invalid value for n_field')
     return
  end select
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Rotation model
  !
  select case (rotation_model)

  case(1)
     call cgyro_info('Rotation terms: O(mach) only (traditional GYRO)')
  case(2)
     call cgyro_info('Rotation terms: O(mach + mach^2) (full Sugama rotation)')
  case(3)
     call cgyro_info('Rotation terms: O(mach^2) only (isolated centrifugal for testing)')
  case(4)
     call cgyro_info('Rotation terms: O(mach^2) GKW CF TRAP term only (for testing)') 
  case(5)
     call cgyro_info('Rotation terms: O(mach^2) GKW CF DRIFT term only (for testing)')   

  case default
     call cgyro_error('Invalid value for rotation_model')
     return

  end select

  !------------------------------------------------------------------------
  ! For adiabatic electrons, n=0 eqn does not have proper electron response
  ! for sonic rotation and may give rise to instability at large Mach
  ! Even for rotation_model=1, the collision eqn does not have the field
  ! correction; explicit_trap_flag=1 will at least do the field correction
  ! for the trapping terms
  if(ae_flag == 1 .and. rotation_model > 1) then
     call cgyro_info('WARNING -- n=0 eqn does not have proper adiabatic ele response')
  endif
  
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Collision model and settings
  !
  ctag(:) = ' '
  ctag(1) = 'x' 

  select case (collision_model)  

  case (1) 
     call cgyro_info('Collision model: Lorentz ee+ei')
  case (2) 
     call cgyro_info('Collision model: Connor')
  case (4) 
     call cgyro_info('Collision model: Sugama')
     ctag(2) = 'x'
  case(5)
     call cgyro_info('Collision model: Simple Lorentz ee+ei')
  case(6)
     call cgyro_info('Collision model: Landau')
  case(7)
     call cgyro_info('Collision model: New (Galerkin) Sugama')
  case default
     call cgyro_error('Invalid value for collision_model')
     return

  end select

  if (collision_model == 6 .or. collision_model == 7) then
     ! if any is below -1 the l-numbers are unlimited.
     if (collision_field_max_l>=-1) then
        write(outstr,*) collision_field_max_l
        call cgyro_info('Field particle collisions limited to l<='//trim(outstr)//' collisions')
     endif
     if (collision_test_max_l>=-1) then
        write(outstr,*) collision_test_max_l
        call cgyro_info('Test particle collisions limited to l<='//trim(outstr)//' collisions')
     endif
  end if
  if (collision_model /= 5) then

     if(collision_model /= 1) then
        select case (collision_mom_restore)
        case(0)
           ! Collision momentum restoring: off       
        case(1)
           ! Collision momentum restoring: on
           ctag(3)= 'x'
        case default
           call cgyro_error('Invalid value for collision_mom_restore')
           return
        end select
     endif

     select case (collision_field_model)
     case(0)
        ! Collision field corrections : off
     case (1)
        ! Collision field corrections : on
        ctag(7) = 'x'
     case default
        call cgyro_error('Invalid value for collision_field_model')
        return
     end select
  endif

  if (collision_model == 4) then

     select case (collision_ene_diffusion)
     case(0)
        ! Collision energy diffusion  : off
     case(1)
        ! Collision energy diffusion  : on
        ctag(2) = 'x'
     case default
        call cgyro_error('Invalid value for collision_ene_diffusion')
        return
     end select

     select case (collision_ene_restore)
     case(0)
        ! Collision energy restoring  : off
     case(1)
        ! Collision energy restoring  : on
        ctag(4) = 'x'
     case default
        call cgyro_error('Invalid value for collision_ene_restore')
        return
     end select

     select case (collision_kperp)
     case(0)
        ! Collision kperp corrections : off
     case(1)
        ! Collision kperp corrections : on
        ctag(5) = 'x'
     case default
        call cgyro_error('Invalid value for collision_kperp')
        return
     end select

  endif

  if (collision_model /= 5 .and. collision_model /= 1) then
     select case (collision_ion_model)
     case(0)
        ! Collisional ions: on
        ctag(6) = 'x'
     case(1)
        ! Collisional ions: off
     case default
        call cgyro_error('Invalid value for collision_ion_model')
        return
     end select
  endif

  if (collision_model > 5) ctag = 'x'
  
  call cgyro_info('Collision terms: L D Rm Re kp ions field')
  call cgyro_info('               '// &
       '  '//ctag(1)// &
       ' '//ctag(2)// &
       ' '//ctag(3)// &
       '  '//ctag(4)// &
       '  '//ctag(5)// &
       '   '//ctag(6)// &
       '     '//ctag(7))
  
  if (collision_model == 5 .or. collision_model == 1) then
     select case(z_eff_method)
     case(1)
        call cgyro_info('Collision model Z_eff: Using Z_eff input')
     case(2)
        call cgyro_info('Collision model Z_eff: Computed from ni and Zi')
     case default
        call cgyro_error('Invalid value for z_eff_method')
        return  
     end select
  endif
  !------------------------------------------------------------------------


  if (rotation_model > 1 .and. collision_model == 5) then
     call cgyro_error('Simple collisions not available with rotation_model > 1')
     return
  endif


  if (collision_test_mode>0) then
     write(unit=outstr,fmt='(A,I1)') 'Warning, 0<collision_test_mode=',collision_test_mode
     call cgyro_info(trim(outstr))
  end if
     
  
  !------------------------------------------------------------------------
  ! Check profile parameters
  !
  do is=1,n_species
     if (dens(is) <= 0.0) then
        call cgyro_error('Densities must be positive.')
        return
     endif
     if (temp(is) <= 0.0) then
        call cgyro_error('Temperatures must be positive.')
        return
     endif
     if (nu(is) < 0.0) then
        call cgyro_error('Collision frequencies must be non-negative.')
        return
     endif
     if (abs(z(is)) < epsilon(0.0)) then
        call cgyro_error('Charges must be non-zero.')
     endif
  enddo
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Dissipation checks
  !
  if (nup_theta < 1 .or. nup_theta > 3) then
     call cgyro_error('Invalid value for nup_theta')
     return
  endif
  if (nup_radial < 1 .or. nup_radial > 4) then
     call cgyro_error('Invalid value for nup_radial')
     return
  endif
  if (nup_alpha < 1 .or. nup_alpha > 4) then
     call cgyro_error('Invalid value for nup_alpha')
     return
  endif
  !------------------------------------------------------------------------

  if (global_flag == 1) then 
     call cgyro_info('##################### IMPORTANT ######################')
     call cgyro_info('#       GLOBAL_FLAG=1 not ready for production       #')
     call cgyro_info('#  See https://github.com/gafusion/gacode/issues/451 #')
     call cgyro_info('######################################################')
  endif

end subroutine cgyro_check
