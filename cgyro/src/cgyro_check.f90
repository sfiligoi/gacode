subroutine cgyro_check

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is
  logical :: lfe

  !-----------------------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_xi,2) /= 0) then 
     call cgyro_error('n_xi must be even.')
     return
  endif

  if (zf_test_flag == 0 .and. modulo(n_radial,box_size) /= 0) then 
     call cgyro_info('RESOLUTION WARNING -- n_radial not a multiple of box_size.')
  endif

  if (n_radial < (n_toroidal-1)*box_size) then
     call cgyro_info('RESOLUTION WARNING -- n_radial < n*box_size')
  endif

  if (n_species > 6) then
     call cgyro_error('n_species <= 6.')
     return
  endif 
 !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Time integration
  !
  if (implicit_flag == 0) then
     call cgyro_info('Integration: RK4 + Implicit C')
  else
     call cgyro_info('Integration: Implicit streaming + RK4 + Implicit C')
  endif
  !------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Profile checks
  !
  select case(profile_model)

  case (1)
     call cgyro_info('Profile model 1: local input (input.cgyro)')

  case (2)
     call cgyro_info('Profile model 2: experimental (input.profiles)')

  case default
     call cgyro_error('Invalid value for profile_model')
     return

  end select
  !-----------------------------------------------------------------------

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
     call cgyro_error('Invalid value for n_field.')
     return
  end select

  if (collision_model == 1 .and. ae_flag == 1) then
     call cgyro_error('Collision_model=1 requires kinetic electrons')
     return
  endif
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Collision model and settings
  !
  select case (collision_model)  

  case (1) 
     call cgyro_info('Collision model 1: Lorentz ee+ei')

  case (2) 
     call cgyro_info('Collision model 2: Connor')

  case (3) 
     call cgyro_info('Collision model 3: Reduced Hirshman-Sigmar')

  case (4) 
     call cgyro_info('Collision model 4: Sugama')

  case default
     call cgyro_error('Invalid value for collision_model')
     return

  end select

  select case (collision_mom_restore)
  case(0)
     call cgyro_info('Collision momentum restoring: off')
  case(1)
     call cgyro_info('Collision momentum restoring: on')
  case default
     call cgyro_error('Invalid value for collision_mom_restore')
     return
  end select

  if (collision_model == 4) then
     select case (collision_ene_diffusion)
     case(0)
        call cgyro_info('Collision energy diffusion  : off')
     case(1)
        call cgyro_info('Collision energy diffusion  : on')
     case default
        call cgyro_error('Invalid value for collision_ene_diffusion')
        return
     end select
     select case (collision_ene_restore)
     case(0)
        call cgyro_info('Collision energy restoring  : off')
     case(1)
        call cgyro_info('Collision energy restoring  : on')
     case default
        call cgyro_error('Invalid value for collision_ene_restore')
        return
     end select
     select case (collision_kperp)
     case(0)
        call cgyro_info('Collision kperp corrections : off')
     case(1)
        call cgyro_info('Collision kperp corrections : on')
     case default
        call cgyro_error('Invalid value for collision_kperp')
        return
     end select
  endif

  select case (collision_field_model)
  case(0)
     call cgyro_info('Collision field corrections : off')
  case (1)
     call cgyro_info('Collision field corrections : on')
  case default
     call cgyro_error('Invalid value for collision_kperp')
     return
  end select
  !
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  

  case (1) 
     call cgyro_info('Equilibrium model 1: s-alpha')
     if (profile_model == 2) then
        call cgyro_error('s-alpha equilibrium model not valid with experimental profiles')
        return
     endif

  case (2) 
     call cgyro_info('Equilibrium model 2: Miller')

  case (3) 
     call cgyro_info('Equilibrium model 3: General (Fourier)')

     if (geo_ny <= 0) then
        call cgyro_error('Fourier geometry coefficients missing.')
        return
     endif

  case default

     call cgyro_error('Invalid value for equilibrium_model')
     return

  end select
  !------------------------------------------------------------------------

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
     if (z(is) == 0) then
        call cgyro_error('Charge must be non-zero.')
        return
     endif
  enddo
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Theta dissipation
  !
  select case (nup_theta)  

  case (1) 
     call cgyro_info('Theta dissipation : 2nd order')

  case (2) 
     call cgyro_info('Theta dissipation : 4th order')

  case (3) 
     call cgyro_info('Theta dissipation : 6th order')

  case default
     call cgyro_error('Invalid value for nup_theta')
     return

  end select
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Radial dissipation
  !
  select case (nup_radial)  

  case (1) 
     call cgyro_info('Radial dissipation: 2nd order')

  case (2) 
     call cgyro_info('Radial dissipation: 4th order')

  case (3) 
     call cgyro_info('Radial dissipation: 6th order')

  case (4) 
     call cgyro_info('Radial dissipation: 8th order')

  case default
     call cgyro_error('Invalid value for nup_radial')
     return

  end select
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Check for existence of restart file
  !
  if (restart_mode == 1) then

     inquire(file=trim(path)//runfile_restart,exist=lfe)
     if (lfe .eqv. .false.) then
        call cgyro_error('Missing restart file')
        return
     endif

  endif
  !------------------------------------------------------------------------

end subroutine cgyro_check
