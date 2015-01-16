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
     call cgyro_error('n_radial must be a multiple of m_box.')
     return
  endif

  if (n_species > 6) then
     call cgyro_error('n_species <= 6.')
     return
  endif
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
     call cgyro_info('Transverse electromagnetic fluctuations (Phi,A_par)')

  case (3)
     call cgyro_info('Compressional electromagnetic fluctuations (Phi,A_par,B_par)')
     stop
  case default
     call cgyro_error('Invalid value for n_field.')
     return
  end select

  if (collision_model == 1 .and. ae_flag == 1) then
     call cgyro_error('collision_model=1 requires kinetic electrons')
     return
  endif
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Collision model

  select case (collision_model)  

  case(0)
     call cgyro_info('Collision model = none')

  case (1) 
     call cgyro_info('Collision model = CONNOR EE+EI LORENTZ only')

  case (2) 
     call cgyro_info('Collision model = Connor')

  case (3) 
     call cgyro_info('Collision model = Reduced Hirshman-Sigmar')

  case (4) 
     call cgyro_info('Collision model = Ad hoc Fokker-Planck')

  case default

     call cgyro_error('Invalid value for collision_model')
     return

  end select

  select case (collision_mom_restore)
  case(0)
     call cgyro_info('Collision momentum restoring = not included')
  case(1)
     call cgyro_info('Collision momentum restoring = included')
  case default
     call cgyro_error('Invalid value for collision_mom_restore')
     return
  end select

  if(collision_model == 4) then
     select case (collision_ene_diffusion)
     case(0)
        call cgyro_info('Collision energy diffusion = not included')
     case(1)
        call cgyro_info('Collision energy diffusion = included')
     case default
        call cgyro_error('Invalid value for collision_ene_diffusion')
        return
     end select
     select case (collision_ene_restore)
     case(0)
        call cgyro_info('Collision energy restoring = not included')
     case(1)
        call cgyro_info('Collision energy restoring = included')
     case default
        call cgyro_error('Invalid value for collision_ene_restore')
        return
     end select
     select case (collision_kperp)
     case(0)
        call cgyro_info('Collision kperp corrections = not included')
     case(1)
        call cgyro_info('Collision kperp corrections = included')
     case default
        call cgyro_error('Invalid value for collision_kperp')
        return
     end select
  endif
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  

  case (0) 
     call cgyro_info('Equlibrium model = s-alpha')

  case (1) 
     call cgyro_error('Invalid value for equilibrium_model.')
     return

  case (2) 
     call cgyro_info('Equlibrium model = Miller')

  case (3) 
     call cgyro_info('Equlibrium model = General (Fourier)')

     if (geo_ny <= 0) then
        call cgyro_error('Fourier geometry coefficients missing.')
        return
     endif

  case default

     call cgyro_error('ERROR: (CGYRO) equilibrium_model invalid')
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
  ! Check for existence of restart file
  !
  if (restart_mode == 1) then

     inquire(file=trim(path)//runfile_restart,exist=lfe)
     if (lfe .eqv. .false.) then
        call cgyro_error('ERROR: (CGYRO) Missing restart file')
        return
     endif

  endif
  !------------------------------------------------------------------------

end subroutine cgyro_check
