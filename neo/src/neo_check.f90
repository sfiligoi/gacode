subroutine neo_check
  
  use neo_globals
  
  implicit none
  
  integer :: ir, is
  
  !-----------------------------------------------------------
  ! Grid parameter checks
  !
  if (n_radial == 1 .and. n_order > 1) then
     print *, 'n_radial=1 must be run with n_order = 1'
     stop
  endif
  !
  if (modulo(n_theta,2) == 0) then 
     print *,'n_theta must be odd.'
     stop
  endif
  !
  if(n_species > 6) then
     print *, 'max n_species is 6'
     stop
  endif
  
  if(rho_in < 0) then
     print *, 'rho_unit must be positive'
     stop
  endif

  !-----------------------------------------------------------
  if(write_out_mode > 1) then
     print *,'SWITCHES'
     print *,'---------------'
  end if
  !------------------------------------------------------------

  ! Write mode
  select case(write_out_mode)
     case(0)
        ! Do not write files, do not print stdout
     case(1)
        ! Do write files, do not print stdout
     case(2)
        ! Do write files, do print stdout
     case default
        print *,'INVALID: write_out_mode'
        stop
  end select

  ! Simulation model
  select case (sim_model)
  case(0)
     if(write_out_mode > 1) then
        print *,'sim_model    : THEORY'
     endif
  case(1)
     if(write_out_mode > 1) then
        print *,'sim_model    : NUMERICAL'
     endif
  case default   
     print *,'INVALID: sim_model'
     stop
  end select

  ! Collision model
  !
  select case (collision_model)  
     
  case (1) 
     
     if(write_out_mode > 1) then
        print *,'collision_model    : CONNOR'
     endif
     
  case (2) 
     
     if(write_out_mode > 1) then
        print *,'collision_model    : REDUCED HIRSHMAN-SIGMAR'
     end if
     
  case (3) 
     
     if(write_out_mode > 1) then
        print *,'collision_model    : FULL HIRSHMAN-SIGMAR'
     endif
     
  case (4) 
     
     if(write_out_mode > 1) then
        print *,'collision_model    : FULL HIRSHMAN-SIGMAR WITH KROOK-LIKE ENERGY DIFFUSION' 
     end if
     
  case default
     
     print *,'INVALID: collision_model'
     stop
     
  end select

  if(collision_model == 3 .or. collision_model == 4) then
     do ir=1, n_radial
        do is=2, n_species
           if(abs(temp(is,ir)-temp(1,ir)) > 1e-3) then
              if(write_out_mode > 1) then
                 print *, 'WARNING: USE OF FULL HIRSHMAN-SIGMAR OPERATOR WITH UNEQUAL TEMPS'
              endif
              exit
           endif
        enddo
     enddo
  endif

  !------------------------------------------------------------
  
  !------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  
     

  case (0) 
     
     if(write_out_mode > 1) then
        print *,'equilibrium_model  : S-ALPHA'
     end if

  case (1) 
     
     if(write_out_mode > 1) then
        print *,'equilibrium_model  : LARGE-ASPECT-RATIO'
     end if
     
  case (2) 
     
     if(write_out_mode > 1) then
        print *,'equilibrium_model  : MILLER'
     end if

  case (3) 
     
     if(write_out_mode > 1) then
        print *,'equilibrium_model  : GENERAL'
     end if   

     if(geo_ny <= 0) then
        print *, 'Geometry coefficients missing'
        stop
     endif

  case default
     
     print *,'INVALID: equilibrium_model'
     stop
     
  end select
  !------------------------------------------------------------
  
  !------------------------------------------------------------
  ! Profile model
  !
  select case (profile_model)  
     
  case (1) 
     
     if(n_radial > 1) then
        print *, 'profile_model=1 must be run with n_radial = 1'
        stop
     endif
     ir=1
     do is=1,n_species
        if(dens(is,ir) <= 0.0) then
           print *, 'density must be positive'
           stop
        end if
        if(temp(is,ir) <= 0.0) then
           print *, 'temperature must be positive'
           stop
        end if
        if(nu(is,ir) <= 0.0) then
           print *, 'collision frequency must be positive'
           stop
        end if
        if(z(is) == 0.0) then
           print *, 'charge must be non-zero'
           stop
        end if
     enddo

     if(write_out_mode > 1) then
        print *,'profile_model      : LOCAL'
     end if
     
  case (2) 
     
     if(write_out_mode > 1) then
        print *,'profile_model      : GLOBAL PROFILE'
     end if
     
     select case (profile_erad0_model)
     case(0)
        if(write_out_mode > 1) then
           print *, 'GLOBAL PROFILE profile_erad0_model: ERAD0 NOT INCLUDED'
        end if
     case (1)
        if(write_out_mode > 1) then
           print *, 'GLOBAL PROFILE profile_erad0_model: ERAD0 INCLUDED'
        end if
     case default
        print *,'INVALID: profile_erad0_model'
        stop
     end select
     
     select case (profile_temprescale_model)
     case(0)
        if(write_out_mode > 1) then
           print *, 'GLOBAL PROFILE profile_temprescale_model:  PROFILE TEMPERATURES ARE NOT RE-SCALED'
        end if
     case (1)
        if(write_out_mode > 1) then
           print *, 'GLOBAL PROFILE profile_temprescale_model: PROFILE TEMPERATURES ARE RE-SCALED TO THE ELECTRON TEMP'
        end if
     case default
        print *,'INVALID: profile_temprescale_model'
        stop
     end select
     
     select case (profile_equilibrium_model)
     case(0)
        if(write_out_mode > 1) then
           print *, 'GLOBAL PROFILE profile_equilibrium_model: WITH S-ALPHA GEOMETRY'
        end if
     case (1)
        if(write_out_mode > 1) then
           print *, 'GLOBAL PROFILE profile_equilibrium_model: WITH MILLER GEOMETRY'
           if(abs(profile_delta_scale-1.0) > epsilon(0.) ) then
              print *, 'GLOBAL PROFILE profile_equilibrium_model: DELTA AND S_DELTA ARE RE-SCALED'
           endif
           if(abs(profile_zeta_scale-1.0) > epsilon(0.) ) then
              print *, 'GLOBAL PROFILE profile_equilibrium_model: ZETA AND S_ZETA ARE RE-SCALED'
           endif
           if(abs(profile_zmag_scale-1.0) > epsilon(0.) ) then
              print *, 'GLOBAL PROFILE profile_equilibrium_model: ZMAG AND S_MAG ARE RE-SCALED'
           endif
        end if
     case (2)
        if(write_out_mode > 1) then
           print *, 'GLOBAL PROFILE profile_equilibrium_model: WITH GENERAL GEOMETRY'
        endif
     case default
        print *,'INVALID: profile_equilibrium_model'
        stop
     end select
     
  case(3)
     if(write_out_mode > 1) then
        print *,'profile_model      : GLOBAL PROFILE TEST'
     end if
     
  case default
     
     print *,'INVALID: profile_model'
     stop
     
  end select

  !-----------------------------------------------------------
  ! Sign of B checks
  if(sign_q > 0.0) then
     if(write_out_mode > 1) then
        print *,'sign_q: POSITIVE'
     end if
  else
     if(write_out_mode > 1) then
        print *, 'sign_q: NEGATIVE'
     end if
  end if
  
  if(sign_bunit > 0.0) then
     if(write_out_mode > 1) then
        print *,'sign_bunit: POSITIVE (BT CW)'
     end if
  else
     if(write_out_mode > 1) then
        print *, 'sign_bunit: NEGATIVE (BT CCW)'
     end if
  end if
  

  !-----------------------------------------------------------
  ! Rotation parameter checks
  !
  select case (rotation_model)     
  case (1)
     if(write_out_mode > 1) then
        print *,'rotation model: ROTATION EFFECTS NOT INCLUDED'
     end if
  case (2)
     if(write_out_mode > 1) then
        print *,'rotation model: ROTATION EFFECTS INCLUDED'
     end if
     if (n_order > 1) then
        print *, 'rotation_model=2 must be run with n_order = 1'
        stop
     endif
     
  case default
     print *,'INVALID: rotation_model'
     stop
  end select
  
  select case (zf_model)
  case(0)
     if(write_out_mode > 1) then
        print *, 'zf_model: NEOCLASSICAL TRANSPORT TEST CASE'
     end if
  case (1)
     if(write_out_mode > 1) then
        print *, 'ZF_model: ZF DAMPING TEST CASE'
     end if
     if(rotation_model == 2) then
        print *, 'zf_model = 1 must be run without rotation effects'
        stop
     endif
     if (n_order > 1) then
        print *, 'zf_model = 1 must be run with n_order = 1'
        stop
     endif
  case default
     print *,'INVALID: zf_model'
     stop
  end select
     
  !------------------------------------------------------------
  
  if(write_out_mode > 1) then
     print *
     print *,'GRID DIMENSIONS'
     print *,'---------------'
     print 10,'n_radial',n_radial
     print 10,'n_energy',n_energy
     print 10,'n_xi',n_xi
     print 10,'n_theta',n_theta
     
     do ir=1,n_radial
        print *
        print *,'PHYSICS PARAMETERS'
        print *,'------------------'
        print 20,'r/R',r(ir)
        print 20,'dphi0/dr',dphi0dr(ir)
        print 20,'omega_rot',omega_rot(ir)
        print 20,'omega_rot_deriv',omega_rot_deriv(ir)
        print 20,'q',q(ir)
        print 20,'s',shat(ir)
        print 20,'shift',shift(ir)
        print 20,'kappa',kappa(ir)
        print 20,'s_kappa',s_kappa(ir)
        print 20,'delta',delta(ir)
        print 20,'s_delta',s_delta(ir)
        print 20,'zeta',zeta(ir)
        print 20,'s_zeta',s_zeta(ir)
        print 20,'zmag',zmag(ir)
        print 20,'s_zmag',s_zmag(ir)

        do is=1,n_species
           print *
           print '(t2,a,i1)','Species ',is
           print *, '----------'
           print 10,'Z',z(is)
           print 20,'dens',dens(is,ir)
           print 20,'temp',temp(is,ir)
           print 20,'mass',mass(is)
           print 20,'rho',rho(ir)*sqrt(temp(is,ir)*mass(is)/temp(1,ir)/mass(1))
           print 20,'a/Ln',dlnndr(is,ir)
           print 20,'a/LT',dlntdr(is,ir)
           print 20,'nu',nu(is,ir)
        enddo
     enddo
 
     print *

  endif
  
10 format(t2,a,':',t14,i2) 
20 format(t2,a,':',t13,1pe12.4) 
  
end subroutine neo_check
