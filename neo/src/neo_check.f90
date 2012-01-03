subroutine neo_check

  use neo_globals

  implicit none

  integer :: ir, is

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
  endif

  !-----------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_theta,2) == 0) then 
     call neo_error('ERROR: (NEO) n_theta must be odd')
     return
  endif
  !
  if(n_species > 6) then
     call neo_error('ERROR: (NEO) max n_species is 6')
     return
  endif

  if(rho_in < 0) then
     call neo_error('ERROR: (NEO) rho_unit must be positive')
     return
  endif

  !-----------------------------------------------------------
  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_neoout,*) 'SWITCHES'
     write(io_neoout,*) '---------------'
  end if
  !------------------------------------------------------------

  ! Simulation model
  select case (sim_model)
  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'sim_model    : THEORY'
     endif
  case(1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'sim_model    : NUMERICAL'
     endif
  case default   
     call neo_error('ERROR: (NEO) invalid sim_model')
     return
  end select

  ! Collision model
  !
  select case (collision_model)  

  case (1) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'collision_model    : CONNOR'
     endif

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'collision_model    : REDUCED HIRSHMAN-SIGMAR'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'collision_model    : FULL HIRSHMAN-SIGMAR'
     endif

  case (4) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'collision_model    : FULL LINEARIZED FOKKER-PLANCK'
     end if

  case (5) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'collision_model    : FULL LINEARIZED FOKKER-PLANCK WITH AD-HOC FIELD PARTICLE TERMS' 
     end if

  case default

     call neo_error('ERROR: (NEO) invalid collision_model')
     return

  end select

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  


  case (0) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'equilibrium_model  : S-ALPHA'
     end if

  case (1) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'equilibrium_model  : LARGE-ASPECT-RATIO'
     end if

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'equilibrium_model  : MILLER'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'equilibrium_model  : GENERAL'
     end if

     if(geo_ny <= 0) then
        call neo_error('ERROR: (NEO) geometry coefficients missing')
        return
     endif

  case default

     call neo_error('ERROR: (NEO) equilibrium_model invalid')
     return

  end select
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Profile model
  !
  select case (profile_model)  

  case (1) 

     if(n_radial > 1) then
        call neo_error('ERROR: (NEO) profile_model=1 must be run with n_radial = 1')
        return
     endif
     ir=1
     do is=1,n_species
        if(dens(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) density must be positive')
           return
        end if
        if(temp(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) temperature must be positive')
           return
        end if
        if(nu(is,ir) <= 0.0) then
           call neo_error('ERROR: (NEO) collision frequency must be positive')
           return
        end if
        if(z(is) == 0.0) then
           call neo_error('ERROR: (NEO) charge must be non-zero')
           return
        end if
     enddo

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'profile_model      : LOCAL'
     end if

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'profile_model      : GLOBAL PROFILE'
     end if

     select case (profile_erad0_model)
     case(0)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_erad0_model: ERAD0 NOT INCLUDED'
        end if
     case (1)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_erad0_model: ERAD0 INCLUDED'
        end if
     case default
        call neo_error('ERROR: (NEO) invalid profile_erad0_model')
        return
     end select

     select case (profile_equilibrium_model)
     case(0)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: WITH S-ALPHA GEOMETRY'
        end if
     case (1)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: WITH MILLER GEOMETRY'
           if(abs(profile_delta_scale-1.0) > epsilon(0.) ) then
              write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: DELTA AND S_DELTA ARE RE-SCALED'
           endif
           if(abs(profile_zeta_scale-1.0) > epsilon(0.) ) then
              write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: ZETA AND S_ZETA ARE RE-SCALED'
           endif
           if(abs(profile_zmag_scale-1.0) > epsilon(0.) ) then
              write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: ZMAG AND S_MAG ARE RE-SCALED'
           endif
        end if
     case (2)
        if(silent_flag == 0 .and. i_proc == 0) then
           write(io_neoout,*) 'GLOBAL PROFILE profile_equilibrium_model: WITH GENERAL GEOMETRY'
        endif
     case default
        call neo_error('ERROR: (NEO) invalid profile_equilibrium_model')
        return
     end select

  case(3)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'profile_model      : GLOBAL PROFILE TEST'
     end if

  case default

     call neo_error('ERROR: (NEO) invalid profile_model')
     return

  end select

  !-----------------------------------------------------------
  ! Sign of B checks
  if(sign_q > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'sign_q: POSITIVE'
     end if
  else
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'sign_q: NEGATIVE'
     end if
  end if

  if(sign_bunit > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'sign_bunit: POSITIVE (BT CW)'
     end if
  else
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'sign_bunit: NEGATIVE (BT CCW)'
     end if
  end if


  !-----------------------------------------------------------
  ! Rotation parameter checks
  !
  select case (rotation_model)     
  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'rotation model: ROTATION EFFECTS NOT INCLUDED'
     end if
  case (2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'rotation model: ROTATION EFFECTS INCLUDED'
     end if
  case default
     call neo_error('ERROR: (NEO) invalid rotation_model')
     return
  end select

  select case (spitzer_model)
  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'spitzer_model: NEOCLASSICAL TRANSPORT TEST CASE'
     end if
  case (1)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_neoout,*) 'spitzer_model: SPITZER TEST CASE'
     end if
  case default
     call neo_error('ERROR: (NEO) invalid spitzer_model')
     return
  end select

  !------------------------------------------------------------

  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_neoout,*)
     write(io_neoout,*) 'GRID DIMENSIONS'
     write(io_neoout,*) '---------------'
     write(io_neoout,10) 'n_radial',n_radial
     write(io_neoout,10) 'n_energy',n_energy
     write(io_neoout,10) 'n_xi',n_xi
     write(io_neoout,10) 'n_theta',n_theta

     do ir=1,n_radial
        write(io_neoout,*) 
        write(io_neoout,*) 'PHYSICS PARAMETERS'
        write(io_neoout,*) '------------------'
        write(io_neoout,20) 'r/R',r(ir)
        write(io_neoout,20) 'dphi0/dr',dphi0dr(ir)
        write(io_neoout,20) 'omega_rot',omega_rot(ir)
        write(io_neoout,20) 'omega_rot_deriv',omega_rot_deriv(ir)
        write(io_neoout,20) 'q',q(ir)
        write(io_neoout,20) 's',shat(ir)
        write(io_neoout,20) 'shift',shift(ir)
        write(io_neoout,20) 'kappa',kappa(ir)
        write(io_neoout,20) 's_kappa',s_kappa(ir)
        write(io_neoout,20) 'delta',delta(ir)
        write(io_neoout,20) 's_delta',s_delta(ir)
        write(io_neoout,20) 'zeta',zeta(ir)
        write(io_neoout,20) 's_zeta',s_zeta(ir)
        write(io_neoout,20) 'zmag',zmag(ir)
        write(io_neoout,20) 's_zmag',s_zmag(ir)

        do is=1,n_species
           write(io_neoout,*) 
           write(io_neoout,'(t2,a,i1)') 'Species ',is
           write(io_neoout,*) '----------'
           write(io_neoout,10) 'Z',z(is)
           write(io_neoout,20) 'dens',dens(is,ir)
           write(io_neoout,20) 'temp',temp(is,ir)
           write(io_neoout,20) 'mass',mass(is)
           write(io_neoout,20) 'a/Ln',dlnndr(is,ir)
           write(io_neoout,20) 'a/LT',dlntdr(is,ir)
           write(io_neoout,20) 'nu',nu(is,ir)
        enddo
     enddo

     write(io_neoout,*)

  endif

  if(silent_flag == 0 .and. i_proc == 0) then
     close(io_neoout)
  endif

10 format(t2,a,':',t14,i2)
20 format(t2,a,':',t13,1pe12.4) 

end subroutine neo_check
