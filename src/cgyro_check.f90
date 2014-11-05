subroutine cgyro_check

  use cgyro_globals

  implicit none

  integer :: is

  if(silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_cgyroout,file=trim(path)//runfile,&
          status='old',position='append')
  endif

  !-----------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_xi,2) /= 0) then 
     call cgyro_error('ERROR: (CGYRO) n_xi must be even')
     return
  endif
  
  if(n_species > 6) then
     call cgyro_error('ERROR: (CGYRO) max n_species is 6')
     return
  endif

  !-----------------------------------------------------------
  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_cgyroout,*) 'SWITCHES'
     write(io_cgyroout,*) '---------------'
  end if
  !------------------------------------------------------------

  ! Collision model
  !
  select case (collision_model)  

  case (0) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'collision_model    : CONNOR EE+EI LORENTZ only'
     endif   

  case (1) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'collision_model    : CONNOR'
     endif

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'collision_model    : REDUCED HIRSHMAN-SIGMAR'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'collision_model    : AD HOC FP'
     endif

  case(-1)
     
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'collision_model    : NONE'
     end if

  case default

     call cgyro_error('ERROR: (CGYRO) invalid collision_model')
     return

  end select

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  


  case (0) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'equilibrium_model  : S-ALPHA'
     end if

  case (1) 

     !if(silent_flag == 0 .and. i_proc == 0) then
     !   write(io_cgyroout,*) 'equilibrium_model  : LARGE-ASPECT-RATIO'
     !end if
     call cgyro_error('ERROR: (CGYRO) equilibrium_model invalid')
     return

  case (2) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'equilibrium_model  : MILLER'
     end if

  case (3) 

     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'equilibrium_model  : GENERAL'
     end if

     if(geo_ny <= 0) then
        call cgyro_error('ERROR: (CGYRO) geometry coefficients missing')
        return
     endif

  case default

     call cgyro_error('ERROR: (CGYRO) equilibrium_model invalid')
     return

  end select
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Check local profile params

  do is=1,n_species
     if(dens(is) <= 0.0) then
        call cgyro_error('ERROR: (CGYRO) density must be positive')
        return
     end if
     if(temp(is) <= 0.0) then
        call cgyro_error('ERROR: (CGYRO) temperature must be positive')
        return
     end if
     if(nu(is) < 0.0) then
        call cgyro_error('ERROR: (CGYRO) collision frequency must be positive')
        return
     end if
     if(z(is) == 0.0) then
        call cgyro_error('ERROR: (CGYRO) charge must be non-zero')
        return
     end if
  enddo

  !------------------------------------------------------------
  ! Toroidal mode number model
  !
  select case (toroidal_model)  

  case(0)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'toroidal_model: Specify k_theta'
     endif
     if(k_theta_rho < 0) then
        call cgyro_error('ERROR: (CGYRO) k_theta must be positive')
        return
     endif
     
  case (1) 
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'toroidal_model: Specify rho'
     endif

  case(2)
     if(silent_flag == 0 .and. i_proc == 0) then
        write(io_cgyroout,*) 'toroidal_model: n=0 test; Specify rho and r_length_rho'
     endif
     
  case default
     
     call cgyro_error('ERROR: (CGYRO) invalid toroidal_model')
     return
     
  end select

  !------------------------------------------------------------

  if(silent_flag == 0 .and. i_proc == 0) then
     write(io_cgyroout,*)
     write(io_cgyroout,*) 'GRID DIMENSIONS'
     write(io_cgyroout,*) '---------------'
     write(io_cgyroout,10) 'n_radial',n_radial
     write(io_cgyroout,10) 'n_energy',n_energy
     write(io_cgyroout,10) 'n_xi',n_xi
     write(io_cgyroout,10) 'n_theta',n_theta
     
     write(io_cgyroout,*) 
     write(io_cgyroout,*) 'PHYSICS PARAMETERS'
     write(io_cgyroout,*) '------------------'
     write(io_cgyroout,*) 'toroidal_num', toroidal_num
     write(io_cgyroout,20) 'r/R',rmin
     write(io_cgyroout,20) 'q',q
     write(io_cgyroout,20) 's',shat
     write(io_cgyroout,20) 'shift',shift
     write(io_cgyroout,20) 'kappa',kappa
     write(io_cgyroout,20) 's_kappa',s_kappa
     write(io_cgyroout,20) 'delta',delta
     write(io_cgyroout,20) 's_delta',s_delta
     write(io_cgyroout,20) 'zeta',zeta
     write(io_cgyroout,20) 's_zeta',s_zeta
     write(io_cgyroout,20) 'zmag',zmag
     write(io_cgyroout,20) 's_zmag',s_zmag

     do is=1,n_species
        write(io_cgyroout,*) 
        write(io_cgyroout,'(t2,a,i1)') 'Species ',is
        write(io_cgyroout,*) '----------'
        write(io_cgyroout,10) 'Z',z(is)
        write(io_cgyroout,20) 'dens',dens(is)
        write(io_cgyroout,20) 'temp',temp(is)
        write(io_cgyroout,20) 'mass',mass(is)
        write(io_cgyroout,20) 'a/Ln',dlnndr(is)
        write(io_cgyroout,20) 'a/LT',dlntdr(is)
        write(io_cgyroout,20) 'nu',nu(is)
     enddo

     write(io_cgyroout,*)

  endif

  if(silent_flag == 0 .and. i_proc == 0) then
     close(io_cgyroout)
  endif

10 format(t2,a,':',t14,i2)
20 format(t2,a,':',t13,1pe12.4) 

end subroutine cgyro_check
