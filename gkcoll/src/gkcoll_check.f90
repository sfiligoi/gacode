subroutine gkcoll_check

  use gkcoll_globals

  implicit none

  integer :: is

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_gkcollout,file=trim(path)//runfile_gkcollout,&
          status='old',position='append')
  endif

  !-----------------------------------------------------------
  ! Grid parameter checks
  !
  !if (modulo(n_theta,2) == 0) then 
  !   call gkcoll_error('ERROR: (GKCOLL) n_theta must be odd')
  !   return
  !endif
  !
  if(n_species > 6) then
     call gkcoll_error('ERROR: (GKCOLL) max n_species is 6')
     return
  endif

  !-----------------------------------------------------------
  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_gkcollout,*) 'SWITCHES'
     write(io_gkcollout,*) '---------------'
  end if
  !------------------------------------------------------------

  ! Collision model
  !
  select case (collision_model)  

  case (1) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'collision_model    : CONNOR'
     endif

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'collision_model    : REDUCED HIRSHMAN-SIGMAR'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'collision_model    : FULL HIRSHMAN-SIGMAR'
     endif

  case (4) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'collision_model    : FULL LINEARIZED FOKKER-PLANCK'
     end if

  case(-1)
     
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'collision_model    : NONE'
     end if

  case (5) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'collision_model    : FULL LINEARIZED FOKKER-PLANCK WITH AD-HOC FIELD PARTICLE TERMS' 
     end if

  case default

     call gkcoll_error('ERROR: (GKCOLL) invalid collision_model')
     return

  end select

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  


  case (0) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'equilibrium_model  : S-ALPHA'
     end if

  case (1) 

     !if(silent_flag == 0 .and. i_proc == 0) then
     !   write(io_gkcollout,*) 'equilibrium_model  : LARGE-ASPECT-RATIO'
     !end if
     call gkcoll_error('ERROR: (GKCOLL) equilibrium_model invalid')
     return

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'equilibrium_model  : MILLER'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'equilibrium_model  : GENERAL'
     end if

     if(geo_ny <= 0) then
        call gkcoll_error('ERROR: (GKCOLL) geometry coefficients missing')
        return
     endif

  case default

     call gkcoll_error('ERROR: (GKCOLL) equilibrium_model invalid')
     return

  end select
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Profile model
  !
  select case (profile_model)  

  case (1) 

     do is=1,n_species
        if(dens(is) <= 0.0) then
           call gkcoll_error('ERROR: (GKCOLL) density must be positive')
           return
        end if
        if(temp(is) <= 0.0) then
           call gkcoll_error('ERROR: (GKCOLL) temperature must be positive')
           return
        end if
        !if(nu(is) <= 0.0) then
        !   call gkcoll_error('ERROR: (GKCOLL) collision frequency must be positive')
        !   return
        !end if
        if(z(is) == 0.0) then
           call gkcoll_error('ERROR: (GKCOLL) charge must be non-zero')
           return
        end if
     enddo

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'profile_model      : LOCAL'
     end if

  case (2) 
     
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'profile_model      : EXPERIMENTAL PROFILE'
     end if
     
  case default
     
     call gkcoll_error('ERROR: (GKCOLL) invalid profile_model')
     return
     
  end select
  
  !-----------------------------------------------------------
  ! Sign of B checks
  if(sign_q > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'sign_q: POSITIVE'
     end if
  else
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'sign_q: NEGATIVE'
     end if
  end if

  if(sign_bunit > 0.0) then
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'sign_bunit: POSITIVE (BT CW)'
     end if
  else
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'sign_bunit: NEGATIVE (BT CCW)'
     end if
  end if

  !------------------------------------------------------------
  ! Toroidal mode number model
  !
  select case (toroidal_model)  

  case(0)
     if(profile_model == 2) then
        call gkcoll_error('ERROR: (GKCOLL) toroidal_model=0 not valid with experimental profiles')
        return
     endif
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'toroidal_model: Specify k_theta'
     endif
     if(k_theta_rho < 0) then
        call gkcoll_error('ERROR: (GKCOLL) k_theta must be positive')
        return
     endif
     
  case (1) 
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'toroidal_model: Specify rho'
     endif
     if(profile_model == 1 .and. rho < 0) then
        call gkcoll_error('ERROR: (GKCOLL) rho_unit must be positive')
        return
     endif

  case(2)
     if(profile_model == 2) then
        call gkcoll_error('ERROR: (GKCOLL) toroidal_model=2 not valid with experimental profiles')
        return
     endif
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_gkcollout,*) 'toroidal_model: n=0 test; Specify rho and r_length_rho'
     endif
     
  case default
     
     call gkcoll_error('ERROR: (GKCOLL) invalid toroidal_model')
     return
     
  end select


  !------------------------------------------------------------

  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_gkcollout,*)
     write(io_gkcollout,*) 'GRID DIMENSIONS'
     write(io_gkcollout,*) '---------------'
     write(io_gkcollout,10) 'n_radial',n_radial
     write(io_gkcollout,10) 'n_energy',n_energy
     write(io_gkcollout,10) 'n_xi',n_xi
     write(io_gkcollout,10) 'n_theta',n_theta
     
     write(io_gkcollout,*) 
     write(io_gkcollout,*) 'PHYSICS PARAMETERS'
     write(io_gkcollout,*) '------------------'
     write(io_gkcollout,*) 'toroidal_num', toroidal_num
     write(io_gkcollout,20) 'r/R',rmin
     write(io_gkcollout,20) 'q',q
     write(io_gkcollout,20) 's',shat
     write(io_gkcollout,20) 'shift',shift
     write(io_gkcollout,20) 'kappa',kappa
     write(io_gkcollout,20) 's_kappa',s_kappa
     write(io_gkcollout,20) 'delta',delta
     write(io_gkcollout,20) 's_delta',s_delta
     write(io_gkcollout,20) 'zeta',zeta
     write(io_gkcollout,20) 's_zeta',s_zeta
     write(io_gkcollout,20) 'zmag',zmag
     write(io_gkcollout,20) 's_zmag',s_zmag

     do is=1,n_species
        write(io_gkcollout,*) 
        write(io_gkcollout,'(t2,a,i1)') 'Species ',is
        write(io_gkcollout,*) '----------'
        write(io_gkcollout,10) 'Z',z(is)
        write(io_gkcollout,20) 'dens',dens(is)
        write(io_gkcollout,20) 'temp',temp(is)
        write(io_gkcollout,20) 'mass',mass(is)
        write(io_gkcollout,20) 'a/Ln',dlnndr(is)
        write(io_gkcollout,20) 'a/LT',dlntdr(is)
        write(io_gkcollout,20) 'nu',nu(is)
     enddo

     write(io_gkcollout,*)

  endif

  if(silent_flag == 0 .and. i_proc == 0) then
     close(io_gkcollout)
  endif

10 format(t2,a,':',t14,i2)
20 format(t2,a,':',t13,1pe12.4) 

end subroutine gkcoll_check
